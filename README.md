<p align="center">
  <img src="PlasLR_logo.png" width="500" title="Final Labelling" alt="Final Labelling">
</p>

# PlasLR: Adaptation of Plasmid Prediction for Error-Prone Long Reads
![GitHub](https://img.shields.io/github/license/anuradhawick/PlasLR)
![GitHub code size in bytes](https://img.shields.io/github/languages/code-size/anuradhawick/PlasLR)

Adaptation of contig based plasmid classification for error prone long reads (PacBio and ONT) to improved genome assembly

## Dependencies
PlasLR is coded purely using C++ (v9) and Python 3.6. To run PlasLR, you will need to install the following python and C++ modules.

### Python dependencies
* numpy 1.16.4 
* scikit-learn 0.22.2.post1
* biopython 1.78
* scipy 1.3.0 
* openTSNE 0.5.1
* umap-learn 0.5.1
* seaborn 0.9.0
* h5py 2.9.0
* tqdm 4.7
* tabulate 0.8

### C++ requirements
* GCC version 9.1.0
* OpenMP 4.5 for multi processing

### Third party programs
* DSK: https://github.com/GATB/dsk
    * Add DSK binaries to your PATH variable

## Downloading PlasLR
To download PlasLR, you have to clone the PlasLR repository to your machine.

```
git clone https://github.com/anuradhawick/PlasLR.git
```

## Compiling the Program

```
cd PlasLR
sh build.sh
```

## Note on Probability Thresholds

PlasLR intends to pick the most accurate classification by either PlasClass or PlasFlow. Users can set the thresholds for selection of plasmids and chromosomes.

* Case: PlasClass

    PlasClass outputs single probability values. Therefore, above `prob_plas` will be plasmids and below `prob_chrom` will be the chromosomes. We recommend leaving these two parameters unspecified so that they can be estimated by PlasLR.

* Case: PlasFlow

    PlasFlow output probabilities under several plasmid and chromosome classes under their phylum. Therefore, the `prob_chrom` and `prob_plas` will be the thresholds for each class. Reads below these thresholds will be treated as unclassified. We recommend leaving these two parameters unspecified so that they can be estimated by PlasLR.

## Running Program

```
usage: PlasLR [-h] --reads READS [--threads THREADS] [--ground-truth <IDS>]
              [--bin-width <bin width for coverage histograms>]
              [--bin-coverage <number of bins for coverage histograms>]
              [--max-mem <Max memory for DSK in Mb>] [--resume] [--plasflow]
              [--plasclass] [--prob-chrom PROB_CHROM] [--prob-plas] [--plots]
              --output OUTPUT

PlasLR Plasmid Classification Corrector

optional arguments:
  -h, --help            show this help message and exit
  --reads READS, -r READS
                        Reads path (FASTQ)
  --threads THREADS, -t THREADS
                        Thread limit
  --ground-truth GROUND_TRUTH, -g GROUND_TRUTH  Read ids of reads (For dry runs with ground truth)
  --bin-width <bin width for coverage histograms>, -bw <bin width for coverage histograms>
                        Value of bw*bs will be the total coverage of k-mers in
                        the coverage histograms. Usually k-mers are shifted
                        towards y-axis due to errors. By defaul b=2; coverages
                        upto 400X
  --bin-coverage BIN_WIDTH, -bc BIN_WIDTH
                        Value of bw*bs will be the total coverage of k-mers in
                        the coverage histograms. Usually k-mers are shifted
                        towards y-axis due to errors. By defaul bc=200;
                        coverages upto 400X
  --max-mem <Max memory for DSK in Mb>, -m <Max memory for DSK in Mb>
                        Default 5000. DSK k-mer counter accepts a max memory
                        parameter. However, the complete pipeline requires
                        5GB+ RAM. This is only to make DSK step faster, should
                        you have more RAM.
  --resume              Continue from the last step or the binning step (which
                        ever comes first). Can save time needed to run DSK and
                        obtain k-mers. Ideal for sensitivity tuning
  --plasflow , -pf      PlasFlow result tsv
  --plasclass , -pc     PlasClass result
  --prob-chrom PROB_CHROM, -C PROB_CHROM
                        Chromosome [Default Auto Detected for plasclass and
                        0.5 for plasflow]
  --prob-plas , -P      Plasmid Threshold [Default Auto Detected for plasclass
                        and 0.7 for plasflow]
  --plots               Whether the initial classifications to be corrected or
                        classify based on already labelled ones
  --output OUTPUT, -o OUTPUT
                        Output directory
```

## Running PlasLR on Test Data

```
cd test_data
tar -xvf data.tar.gz 

python ../PlasLR -r reads.fasta -t 8 --ground-truth truth.txt -pc pc -o plaslr_output --plots

```

## Notes

* We recommend not to use plots argument if you're running on a larger dataset (More than 100000 reads). As it might make the execution take a longer time. 

## Issues

Please submit issues on the github repo. Thanks for using PlasLR.
