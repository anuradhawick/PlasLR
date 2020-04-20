import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
from sklearn.decomposition import PCA
from sklearn.neighbors import KNeighborsClassifier
from sklearn.neighbors import NearestNeighbors
from sklearn.metrics import recall_score
from sklearn.metrics import precision_score
from sklearn.metrics import f1_score
import argparse
import os
from collections import defaultdict
import sys

plt.rcParams.update({'font.size': 16})

palette ={"Plasmid":"C0","Chromosome":"C1","Unclassified":"C2"}

parser = argparse.ArgumentParser(description="""Plasmid Correction""")

parser.add_argument('-pf',
                    metavar='',
                    help="PlasFlow result tsv",
                    type=str,
                    default=None,
                    required=False)
parser.add_argument('-pc',
                    metavar='',
                    help="PlasClass result",
                    type=str,
                    default=None,
                    required=False)
parser.add_argument('-i',
                    metavar='',
                    help="Read ids of reads (For dry runs with ground truth)",
                    type=str,
                    default=None,
                    required=False)
parser.add_argument('-o',
                    metavar='',
                    help="",
                    type=str,
                    required=True)                    
parser.add_argument('-p3',
                    metavar='',
                    help="",
                    type=str,
                    required=True)
parser.add_argument('-p15',
                    metavar='',
                    help="",
                    type=str,
                    required=True)
parser.add_argument('-t',
                    metavar='',
                    help="Threads",
                    type=int,
                    default=8,
                    required=False)

args = parser.parse_args()
classification = args.pf

if classification == None:
    classification = args.pc

if classification == None:
    print("Must provide either -pf or -pc as the input")
    print("Program exitting")
    sys.exit(1)

p3 = args.p3
p15 = args.p15
truth = args.i 
output = args.o
threads = args.t

if not os.path.exists(output):
    os.makedirs(output)

def evalCs(truth, clusters):
    string_output = ""
    clusterIds = {}
    speciesIds = {}
    mat = {}
    t_arr = []
    p_arr = []
    for t, c in zip(truth, clusters):
        species = t
        cluster = c

        t_arr.append(t)
        p_arr.append(c)

        if cluster not in clusterIds:
            clusterIds[cluster] = len(clusterIds)
            for k, v in mat.items():
                v.append(0)
        if species not in mat:
            mat[species] = [0 for x in range(len(clusterIds))]
            speciesIds[species] = len(speciesIds)

        mat[species][clusterIds[cluster]] += 1

    print(" " * (31), "".join(map(lambda x: " " * (20 - len(x)) + x, clusterIds.keys())))
    string_output += " " * (31) + "" + "".join(map(lambda x: " " * (20 - len(x)) + x, clusterIds.keys())) + "\n"
    
    for k, v in mat.items():
        string_output += str(k) + " " + " " * (30 - len(k))  + "".join(map(lambda x: " " * (20 - len(str(x))) + str(x), v)) + "\n"
        print(k, " " * (30 - len(k)),
              "".join(map(lambda x: " " * (20 - len(str(x))) + str(x), v)))


    arrMat = []

    for k, v in mat.items():
        arrMat.append(v)

    arrMatTranspose = [[arrMat[j][i] for j in range(len(arrMat))] for i in range(len(arrMat[0]))]

    rMax = 0
    tot = 0
    for x in arrMat:
        tot += sum(x)
        rMax += max(x)

    precision = precision_score(t_arr, p_arr, average='micro')
    recall = recall_score(t_arr, p_arr, average='micro')
    f1 = f1_score(t_arr, p_arr, average='micro')
    print("\n\n")
    print("Precision ", precision)
    string_output += "\n"
    string_output += "Precision {}\n".format(precision)
    print("Recall ", recall)
    string_output += "Recall    {}\n".format(recall)
    print("F1 Score ", f1)
    string_output += "F1 Score  {}\n".format(f1)
    return string_output

def processLabel(x):
    x = x.replace(">", "")
    if x.lower() in "chromosomes" or "chromosome" in x.lower():
        return "Chromosome"
    elif x.lower() in "plasmids" or "plasmid" in x.lower():
        return "Plasmid"
    else:
        return "Unclassified"

original_labels = []
if args.pf != None:
    # PlasFlow parser
    with open(classification) as clusters:
        clusters = clusters.read().strip().split("\n")
        bar = clusters[0].strip().split("\t")
        _arr = []
        for c in clusters[1:]:
            label = c.strip().split("\t")[5].split(".")[0]
            cs = [float(x) for x in c.strip().split("\t")[6:]]
            am = np.argmax(cs)
            prob = max(cs)
            cluster = label
            
            original_labels.append(label)
            
            cluster = bar[5::][am].split(".")[0]
            cluster = processLabel(cluster)
            
            if prob < 0.7 and cluster.lower() in "plasmid":
                cluster = "Unclassified"
            elif prob < 0.5 and cluster.lower() in "chromosome":
                cluster = "Unclassified"
                
            _arr.append(cluster)
        classification = _arr
# parsing plasclass
elif args.pc != None:
    with open(classification) as clusters:
        clusters = clusters.read().strip().split("\n")
        bar = clusters[0].strip().split("\t")
        _arr = []
        for c in clusters:
            prob = float(c.strip().split("\t")[-1])
    
            if prob > 0.5:
                original_labels.append("Plasmid")
            else:
                original_labels.append("Chromosome")

            if prob > 0.8:
                cluster = "Plasmid"
            elif prob < 0.5:
                cluster = "Chromosome"
            else:
                cluster = "Unclassified"
                
            _arr.append(cluster)
        classification = _arr
if truth:
    truth = [processLabel(x) for x in open(truth).read().split()]
p3 = np.array([[float(y) for y in x.strip().split()] for x in open(p3).read().strip().split("\n")])
p15 = np.array([[float(y) for y in x.strip().split()] for x in open(p15).read().strip().split("\n")])
classification = np.array(classification)
p15_p3 = np.append(p15, p3, axis=1)
pca1 = PCA(n_components=2)
p15_2d = pca1.fit_transform(p15)
pca2 = PCA(n_components=2)
p15_p3_2d = pca2.fit_transform(p15_p3)
pca3 = PCA(n_components=2)
p3_2d = pca3.fit_transform(p3)

# fig = plt.figure(figsize=(5, 5))
# sns.scatterplot(p15_2d[:,0], p15_2d[:,1], hue=truth, palette=palette)
# plt.xlabel("PCA 1")
# plt.ylabel("PCA 2")
# plt.savefig(output + "/truth_cov_comp.png", dpi=1200)

# fig = plt.figure(figsize=(5, 5))
# sns.scatterplot(p15_p3_2d[:,0], p15_p3_2d[:,1], hue=truth, palette=palette)
# plt.xlabel("PCA 1", fontsize=20)
# plt.ylabel("PCA 2", fontsize=20)
# # plt.savefig(output + "/truth_cov_comp.svg", format="svg", bbox_inches='tight')
# plt.savefig(output + "/truth_cov_comp.eps", dpi=1200, bbox_inches='tight')
# plt.savefig(output + "/truth_cov_comp.png", dpi=300, bbox_inches='tight')

# fig = plt.figure(figsize=(5, 5))
# sns.scatterplot(p15_p3_2d[:,0], p15_p3_2d[:,1], hue=classification, palette=palette)
# plt.xlabel("PCA 1", fontsize=20)
# plt.ylabel("PCA 2", fontsize=20)
# plt.savefig(output + "/labels_before_removal.svg", format="svg", bbox_inches='tight')
# plt.savefig(output + "/labels_before_removal.eps", dpi=1200, bbox_inches='tight')
# plt.savefig(output + "/labels_before_removal.png", dpi=300, bbox_inches='tight')

# preliminary label removal of non-confident ones
nbrs = NearestNeighbors(n_neighbors=50, algorithm='ball_tree', n_jobs=threads).fit(p15_p3_2d)
distances, indices = nbrs.kneighbors(p15_p3_2d)

classification2 = list(classification)

for i in indices:
    _label = classification[i[0]]
    _others = list(map(lambda x: classification[x], i[1:]))
    _vote = 0
    _all_votes = 0
    for x in _others:
        if x == "Unclassified":
            continue
        if x == _label:
            _vote += 1
        _all_votes += 1
    
    if _vote/max(_all_votes, 1) < 0.1 and _all_votes > 30:
        classification2[i[0]] = "Unclassified"
    # if _vote/max(_all_votes, 1) < 0.2 and _all_votes > 30 and _label == "Chromosome":
    #     classification2[i[0]] = "Unclassified"

# fig = plt.figure(figsize=(5, 5))
# sns.scatterplot(p15_p3_2d[:,0], p15_p3_2d[:,1], hue=classification2, palette=palette)
# plt.xlabel("PCA 1", fontsize=20)
# plt.ylabel("PCA 2", fontsize=20)
# plt.savefig(output + "/labels_after_removal.svg", format="svg", bbox_inches='tight')
# plt.savefig(output + "/labels_after_removal.eps", dpi=1200, bbox_inches='tight')
# plt.savefig(output + "/labels_after_removal.png", dpi=300, bbox_inches='tight')


# separation of labelled and unlabelled nodes
labelled_p15_2d = []
labelled_p3_2d = []
labelled_p15_p3_2d = []
labels_p15 = []
unlabelled_p15_2d = []
unlabelled_p3_2d = []
unlabelled_p15_p3_2d = []
truth_unlb = []
truth_lbl = []
read_tracker = [x for x in range(len(p3_2d))]
r_id_labelled = []
r_id_unlabelled = []

if truth:
    for _label, _p15, _p3, _p15_p3, _t, _rid in zip(classification2, p15_2d, p3_2d, p15_p3_2d, truth, read_tracker):
        if _label == "Unclassified":
            unlabelled_p15_2d.append(_p15)
            unlabelled_p3_2d.append(_p3)
            unlabelled_p15_p3_2d.append(_p15_p3)
            truth_unlb.append(_t)
            r_id_unlabelled.append(_rid)
        else:
            labels_p15.append(_label)
            labelled_p15_2d.append(_p15)
            labelled_p3_2d.append(_p3)
            labelled_p15_p3_2d.append(_p15_p3)
            truth_lbl.append(_t)
            r_id_labelled.append(_rid)
else:
    for _label, _p15, _p3, _p15_p3, _rid in zip(classification2, p15_2d, p3_2d, p15_p3_2d, read_tracker):
        if _label == "Unclassified":
            unlabelled_p15_2d.append(_p15)
            unlabelled_p3_2d.append(_p3)
            unlabelled_p15_p3_2d.append(_p15_p3)
            r_id_unlabelled.append(_rid)
        else:
            labels_p15.append(_label)
            labelled_p15_2d.append(_p15)
            labelled_p3_2d.append(_p3)
            labelled_p15_p3_2d.append(_p15_p3)
            r_id_labelled.append(_rid)

labelled_p15_2d = np.array(labelled_p15_2d)
labelled_p3_2d = np.array(labelled_p3_2d)
labelled_p15_p3_2d = np.array(labelled_p15_p3_2d)
labels_p15 = np.array(labels_p15)
unlabelled_p15_2d = np.array(unlabelled_p15_2d)
unlabelled_p3_2d = np.array(unlabelled_p3_2d)
unlabelled_p15_p3_2d = np.array(unlabelled_p15_p3_2d)

clf = KNeighborsClassifier(100, weights='uniform', n_jobs=threads)
trained_model = clf.fit(labelled_p15_p3_2d, labels_p15)

imputed_values = trained_model.predict(unlabelled_p15_p3_2d)

all_2d = np.append(labelled_p15_p3_2d, unlabelled_p15_p3_2d, axis=0)
all_labels = np.append(labels_p15, imputed_values, axis=0)
if truth:
    truth_all = np.append(truth_lbl, truth_unlb, axis=0)
r_id_all = r_id_labelled + r_id_unlabelled

# fig = plt.figure(figsize=(5, 5))
# sns.scatterplot(all_2d[:,0], all_2d[:,1], hue=all_labels, palette=palette)
# plt.xlabel("PCA 1", fontsize=20)
# plt.ylabel("PCA 2", fontsize=20)
# plt.savefig(output + "/label_corrected.eps", dpi=1200, bbox_inches='tight')
# plt.savefig(output + "/label_corrected.png", dpi=300, bbox_inches='tight')
# plt.savefig(output + "/label_corrected.svg", format="svg", bbox_inches='tight')

if truth:
    o1 = evalCs(truth, classification)
    o2 = evalCs(truth_all, all_labels)


    with open(output + "/logs.txt", "w+") as f:
        f.write("Performance before correction\n")
        f.write(o1)
        f.write("\n\n")
        f.write("Performance after correction\n")
        f.write(o2)

readIds_classified = list(zip(r_id_all, all_labels))
readIds_classified.sort(key=lambda x: x[0])
with open(output + "/classification.txt", "w+") as f:
    for i, l in readIds_classified:
        f.write("{}\t{}\n".format(i, l))