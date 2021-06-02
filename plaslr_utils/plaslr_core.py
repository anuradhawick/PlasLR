import numpy as np
import argparse
import os
import sys
import logging

from tqdm import tqdm

import seaborn as sns
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

from sklearn.neighbors import KNeighborsClassifier
from sklearn.neighbors import NearestNeighbors
from sklearn.preprocessing import StandardScaler
from sklearn.manifold import TSNE
from sklearn.metrics import precision_score, recall_score
from sklearn.metrics.cluster import adjusted_rand_score

from tabulate import tabulate

logger = logging.getLogger('PlasLR')

palette ={"Plasmid":"C0","Chromosome":"C1", "plasmid":"C0","chromosome":"C1","unclassified":"C2"}

plt.rcParams.update({'font.size': 16})    


def eval_performance(truth, clusters):
    truth = list(map(lambda x: x.lower(), truth))
    clusters = list(map(lambda x: x.lower(), clusters))
    
    plasmid_count = truth.count("plasmid")
    chromosome_count = truth.count("chromosome")
    
    str_output = f"\nPlasmids    = {plasmid_count: 10}\n"
    str_output += f"Chromosomes = {chromosome_count: 10}\n\n"
    
    truth = np.array(truth)
    clusters = np.array(clusters)

    classified_truths = truth[clusters!='unclassified']
    classified_clusters = clusters[clusters!='unclassified']
    c_as_c = 0
    c_as_p = 0
    p_as_p = 0
    p_as_c = 0
    c_as_u = 0
    p_as_u = 0
    
    for t, c in zip(truth, clusters):
        if c==t=='plasmid':
            p_as_p+=1
        elif c==t=='chromosome':
            c_as_c+=1
        elif t.lower()=='chromosome':
            if c.lower()=='plasmid':
                c_as_p+=1
            else:
                c_as_u+=1
        else:
            if c.lower()=='chromosome':
                p_as_c+=1
            else:
                p_as_u+=1
                
    str_output += tabulate([["Chromosome",c_as_c, c_as_p, c_as_u],
                    ["Plasmid",p_as_c, p_as_p, p_as_u]], 
                   headers=["", "Chromosome", "Plasmid", "Unclassified"], tablefmt="fancy_grid")
    str_output += "\n\n"
    pr = precision_score(classified_truths, classified_clusters, average='micro')
    re = recall_score(truth, clusters, average='micro')
    f1 = 2*pr*re/(pr+re)
    str_output += f"Precision   = {(100*pr):3.2f}\n"
    str_output += f"Recall      = {(100*re):3.2f}\n"
    str_output += f"F1 Score    = {(100*f1):3.2f}\n"
    str_output += f"Recovery    = {(100*p_as_p/(p_as_p + p_as_c + p_as_u)):3.2f}\n"
    
    
    return str_output

def get_thresholds(probs):
    total_reads = len(probs)

    for t in np.arange(0.90, 0.7, -0.01):
        t = round(t, 2)
        percentage_plas = sum(np.where(probs > t, 1, 0))/total_reads * 100

        if percentage_plas > 5:
            break

    upper_thresh = t
    for t in np.arange(0.01, 0.5, 0.01):
        t = round(t, 2)
        percentage_chrom = sum(np.where(probs < t, 1, 0))/total_reads * 100

        if percentage_chrom > 20:
            break

    lower_thresh = t

    return lower_thresh, upper_thresh

def tsne_sampled_func(X, threads):
    indices = np.random.permutation(list(range(X.shape[0])))
    reverse = np.argsort(indices)

    x_sample, x_rest = X[indices[:25000]], X[indices[25000:]]
    
    sample_affinities = affinity.PerplexityBasedNN(
        x_sample,
        perplexity=500,
        method="approx",
        n_jobs=8,
        random_state=0,
    )
    
    sample_init = initialization.pca(x_sample, random_state=42)

    sample_embedding = TSNEEmbedding(
        sample_init,
        sample_affinities,
        negative_gradient_method="fft",
        n_jobs=threads,
        verbose=False
    )
    
    sample_embedding1 = sample_embedding.optimize(n_iter=250, exaggeration=12, momentum=0.5)
    sample_embedding2 = sample_embedding1.optimize(n_iter=750, exaggeration=1, momentum=0.8)
    rest_init = sample_embedding2.prepare_partial(x_rest, k=1, perplexity=1/3)
    init_full = np.vstack((sample_embedding2, rest_init))[reverse]
    
    return init_full

def get_best_bins(size_hist, kmers):
    bins = np.zeros(size_hist)
    best_flux = 0
    best_threshold = 0
    best_bins = []
    
    bins_o = np.zeros(size_hist)
    
    for count in tqdm(kmers, desc='Building kmer histogram'):
        if 0 < count < len(bins_o):
            bins_o[count] += 1
    
    for t in range(100, 1000, 10):
        bins = np.where(bins_o > t, 1, 0)
        clean_width = 10
        flux = 0
        
        for i in range(1, size_hist):
            if bins[i-1] != bins[i]:
                flux+=1
        
        if flux > best_flux:
            best_flux = flux
            best_threshold = t
            best_bins = bins
            
    return best_threshold, best_flux

def run_plasmid_correction(p3, p15, readIds, kmer_counts, output, *, threads=8, truth=None, prob_plas=None, prob_chrom=None, plots=False, plasclass=None, plasflow=None, dimension_reduction=None):
    if plots and not os.path.isdir(f"{output}/images/"):
        os.makedirs(f"{output}/images/")

    lines = 0
    logger.info("Loading the data.")
    if plasclass:
        with open(plasclass) as pc_file:
            for line in pc_file:
                if len(line.strip()) > 0:
                    lines += 1

            pc_file.seek(0)
            probs = [float(line.strip().split("\t")[-1]) for line in tqdm(pc_file, total=lines, desc="Loading plasclass results")]
    elif plasflow:
        with open(plasflow) as pf_file:
            lines = -1
            for line in pf_file:
                if len(line.strip()) > 0:
                    lines += 1

            pf_file.seek(0)
            probs = [sum(map(float, line.strip().split("\t")[24:])) for line in tqdm(pf_file.read().strip().split("\n")[1:], total=lines, desc="Loading plasflow results")]   

    probs = np.array(probs)


    if truth is not None:
        truth = np.array([line.strip().replace(">", "") for line in tqdm(open(truth), total=lines, desc="Loading ground truth")])

    p3 = np.array([[float(y) for y in line.strip().split()] for line in tqdm(open(p3), total=lines, desc="Loading composition profiles")])
    p15 = np.array([[float(y) for y in line.strip().split()] for line in tqdm(open(p15), total=lines, desc="Loading coverage profiles")])
    readIds = np.array([line.strip() for line in tqdm(open(readIds), total=lines, desc="Loading read ids")])
    logger.info("Finish loading the data.")


    logger.info("Computing the probability thresholds.")
    if prob_chrom is None and prob_plas is None:
        lower_thresh, upper_thresh = get_thresholds(probs)
    else:
        lower_thresh, upper_thresh = prob_chrom, prob_plas
    logger.debug(f"Lower threshold = {lower_thresh} Upper threshold = {upper_thresh}")
    logger.info("Finish computing the probability thresholds.")

    classification = np.array(list(map(lambda x: "plasmid" if x > upper_thresh else ("chromosome" if x < lower_thresh else "unclassified"), probs)))
    
    if dimension_reduction == 'umap':
        import umap
        mapper = umap.UMAP()
        profiles = mapper.fit_transform(np.concatenate([p15, p3], axis=1))
    elif dimension_reduction == 'tsne':
        from openTSNE import TSNE, TSNEEmbedding, affinity, initialization
        profiles = tsne_sampled_func(np.concatenate([p15, p3], threads))
    else:
        from sklearn.decomposition import PCA
        profiles = PCA(n_components=2).fit_transform(np.concatenate([p15, p3], axis=1))
        
    if plots and truth is not None:
        fig = plt.figure(figsize=(20, 10))
        ax = fig.add_subplot(1, 2, 1)
        sns.scatterplot(profiles[:,0], profiles[:,1], hue=classification, palette=palette, alpha=0.1)
        ax = fig.add_subplot(1, 2, 2)
        sns.scatterplot(profiles[:,0], profiles[:,1], hue=truth, palette=palette, alpha=0.1)
        plt.savefig(f"{output}/images/figure.classification-vs-truth.png", dpi=100, format="png")

    elif plots:
            fig = plt.figure(figsize=(10, 10))
            sns.scatterplot(profiles[:,0], profiles[:,1], hue=classification, palette=palette, alpha=0.1)
            plt.savefig(f"{output}/images/figure.classification.png", dpi=100, format="png")

    # preliminary label removal of non-confident ones
    # nbrs = NearestNeighbors(n_neighbors=50, algorithm='ball_tree', n_jobs=threads).fit(profiles)
    # distances, indices = nbrs.kneighbors(profiles)

    classification_corrected = list(classification)

    # for i in indices:
    #     _label = classification[i[0]]
        
    #     if _label.lower() == "unclassified":
    #         continue
            
    #     _others = list(map(lambda x: classification[x], i[1:]))
    #     _vote = 0
    #     _all_votes = 0
        
    #     for x in _others:
    #         if x.lower() == "unclassified":
    #             continue
    #         if x == _label:
    #             _vote += 1
    #         _all_votes += 1

    #     if _vote/max(_all_votes, 1) < 0.1 and _all_votes > 30:
    #         classification_corrected[i[0]] = "unclassified"
    classification_corrected = np.array(classification_corrected)
    
    if plots:
        fig = plt.figure(figsize=(20, 10))
        ax = fig.add_subplot(1, 2, 1)
        sns.scatterplot(profiles[:,0], profiles[:,1], hue=classification, palette=palette, alpha=0.1)
        ax = fig.add_subplot(1, 2, 2) 
        sns.scatterplot(profiles[:,0], profiles[:,1], hue=classification_corrected, palette=palette, alpha=0.1)

        plt.savefig(f"{output}/images/figure.classification-before-and-after-correction.png", dpi=100, format="png")

    classified_labels = classification_corrected[classification_corrected!="unclassified"]
    classified_profiles = profiles[classification_corrected!="unclassified"]
    classified_readIds = readIds[classification_corrected!="unclassified"]

    unclassified_labels = classification_corrected[classification_corrected=="unclassified"]
    unclassified_profiles = profiles[classification_corrected=="unclassified"]
    unclassified_readIds = readIds[classification_corrected=="unclassified"]

    logger.debug(f"Classified data size = {len(classified_labels)}")
    logger.debug(f"Unclassified data size = {len(unclassified_labels)}")

    knn_clf = KNeighborsClassifier(100, weights='distance', n_jobs=threads)
    trained_model = knn_clf.fit(classified_profiles, classified_labels)
    plasmid_class_idx = list(knn_clf.classes_).index("plasmid")

    predicted_labels = trained_model.predict(unclassified_profiles)

    all_profiles = np.append(classified_profiles, unclassified_profiles, axis=0)
    all_labels = np.append(classified_labels, predicted_labels, axis=0)
    all_probabilities = knn_clf.predict_proba(all_profiles)
    all_plasmid_probabilities = all_probabilities.T[plasmid_class_idx]
    all_read_ids = np.append(classified_readIds, unclassified_readIds, axis=0)

    if truth is not None:
        classified_truth = truth[classification_corrected!="unclassified"]
        unclassified_truth = truth[classification_corrected=="unclassified"]
        all_truth = np.append(classified_truth, unclassified_truth, axis=0)

        logger.info("Performance before PlasLR correction and classification")
        logger.info(eval_performance(truth, classification))
        logger.info("Performance before PlasLR correction and classification")
        logger.info(eval_performance(all_truth, all_labels))
    
    if plots and truth is not None:
        fig = plt.figure(figsize=(20, 10))
        ax = fig.add_subplot(1, 2, 1)
        sns.scatterplot(all_profiles[:,0], all_profiles[:,1], hue=all_labels, palette=palette, alpha=0.1)
        ax = fig.add_subplot(1, 2, 2)
        sns.scatterplot(all_profiles[:,0], all_profiles[:,1], hue=all_truth, palette=palette, alpha=0.1)
        plt.savefig(f"{output}/images/figure.final-vs-truth.png", dpi=100, format="png")
    elif plots:
        fig = plt.figure(figsize=(10, 10))
        sns.scatterplot(all_profiles[:,0], all_profiles[:,1], hue=all_labels, palette=palette, alpha=0.1)
        plt.savefig(f"{output}/images/figure.final.png", dpi=100, format="png")

    logger.info(f"Writing results to {output}/final.txt")
    with open(f"{output}/final.txt", "w+") as final_result:
        for readId, label, probability in tqdm(zip(all_read_ids, all_labels, all_plasmid_probabilities), desc='Writing results', total=len(all_read_ids)):
            final_result.write(f"{readId}\t{label}\t{probability}\n")
    logger.info(f"Writing results to {output}/final.txt completed")

     
    

if __name__ == '__main__':   
    parser = argparse.ArgumentParser(description="""Plasmid Correction""")

    parser.add_argument('--plasflow', '-pf',
                        help="PlasFlow result tsv",
                        type=str,
                        default=None,
                        required=False)
    parser.add_argument('--plasclass', '-pc',
                        help="PlasClass result",
                        type=str,
                        default=None,
                        required=False)
    parser.add_argument('--ground-truth',
                        help="Read ids of reads (For dry runs with ground truth)",
                        type=str,
                        default=None,
                        required=False)
    parser.add_argument('--read-ids', '-r',
                        help="Read ids of reads (For dry runs with ground truth)",
                        type=str,
                        default=None,
                        required=True)
    parser.add_argument('--output', '-o',
                        help="Output folder path (A new folder with name will be created if not there)",
                        type=str,
                        required=True)                    
    parser.add_argument('--composition-profiles', '-p3',
                        type=str,
                        required=True)
    parser.add_argument('--coverage-profiles', '-p15',
                        type=str,
                        required=True)
    parser.add_argument('--kmer-counts', '-kc',
                        type=str,
                        required=True)
    parser.add_argument('--prob-chrom', '-C',
                        help="Chromosome [Default computed]",
                        type=float,
                        default=None,
                        required=False)
    parser.add_argument('--prob-plas', '-P',
                        help="Plasmid Threshold [Default computed]",
                        type=str,
                        default=None,
                        required=False)
    parser.add_argument('--plots', '-p',
                        action='store_true',
                        help="Whether the initial classifications to be corrected or classify based on already labelled ones")         
    parser.add_argument('--threads', '-t',
                        help="Threads",
                        type=int,
                        default=8,
                        required=False)     

    args = parser.parse_args()
    classification = args.plasflow

    if args.plasclass != None and args.plasflow != None:
        print("Must provide only one of -pf or -pc as the input")
        print("Program exitting")
        sys.exit(1)

    if classification == None:
        classification = args.plasclass

    if classification == None:
        print("Must provide either -pf or -pc as the input")
        print("Program exitting")
        sys.exit(1)

    p3 = args.composition_profiles
    p15 = args.coverage_profiles
    truth = args.ground_truth
    output = args.output
    threads = args.threads
    prob_plas = args.prob_plas
    prob_chrom = args.prob_chrom
    plots = args.plots
    readIds = args.read_ids
    kmer_counts = args.kmer_counts

    if not (prob_plas==None and prob_chrom==None) and (prob_plas!=None or prob_chrom!=None):
        print("In manual mode both prob_plas and prob_chrom must be specified")
        print("Please consider leaving the values for automatic thrshold setting")
        print("Program exitting")
        sys.exit(1)
    
    run_plasmid_correction(p3, p15, readIds, kmer_counts, output, threads=8, truth=truth, prob_plas=prob_plas, prob_chrom=prob_chrom, plots=plots, plasflow=args.plasflow, plasclass=args.plasclass)
