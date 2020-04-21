# PlasLR
Adaptation of contig based plasmid classification for error prone long reads (PacBio and ONT) to improved genome assembly

## Compiling the Program

```
sh build.sh
```

## Running Program

```
usage: PlasLR [-h] -r <READS PATH> [-t <THREADS>] [-i <IDS>]
              [-b <bin width for coverage histograms>]
              [-m <Max memory for DSK in Mb>] [--resume] [-pf] [-pc]
              [-prob_chrom] [-prob_plas] [-label_corr] [-plots] -o <DEST>

PlasLR Plasmid Classification Corrector

optional arguments:
  -h, --help            show this help message and exit
  -r <READS PATH>       Reads path (FASTQ)
  -t <THREADS>          Thread limit
  -i <IDS>              Read ids of reads (For dry runs with ground truth)
  -b <bin width for coverage histograms>
                        Value of bx32 will be the total coverage of k-mers in
                        the coverage histograms. Usually k-mers are shifted
                        towards y-axis due to errors. By defaul b=10;
                        coverages upto 320X
  -m <Max memory for DSK in Mb>
                        Default 5000. DSK k-mer counter accepts a max memory
                        parameter. However, the complete pipeline requires
                        5GB+ RAM. This is only to make DSK step faster, should
                        you have more RAM.
  --resume              Continue from the last step or the binning step (which
                        ever comes first). Can save time needed to run DSK and
                        obtain k-mers. Ideal for sensitivity tuning
  -pf                   PlasFlow result tsv
  -pc                   PlasClass result
  -prob_chrom           Chromosome [Default 0.7 for plasclass and 0.5 for
                        plasflow]
  -prob_plas            Plasmid Threshold [Default 0.3 for plasclass and 0.7
                        for plasflow]
  -label_corr           Whether the initial classifications to be corrected or
                        classify based on already labelled ones
  -plots                Whether the initial classifications to be corrected or
                        classify based on already labelled ones
  -o <DEST>             Output directory
```