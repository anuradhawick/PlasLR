mkdir ./reads

simlord -t 0.1 -mr 1500 -rr ./NC_000906_plasmid.faa --no-sam ./reads/NC_000906_plasmid -c 100
simlord -t 0.1 -mr 1500 -rr ./NC_000914_plasmid.faa --no-sam ./reads/NC_000914_plasmid -c 300
simlord -t 0.1 -mr 1500 -rr ./NC_000918_chromosome.faa --no-sam ./reads/NC_000918_chromosome -c 20
simlord -t 0.1 -mr 1500 -rr ./NC_002570_chromosome.faa --no-sam ./reads/NC_002570_chromosome -c 30

awk '{if(NR%4==1){print ">plasmid"}else if(NR%4==2){print}}' ./reads/NC_000906_plasmid.fastq > ./reads/NC_000906_plasmid.fasta
awk '{if(NR%4==1){print ">plasmid"}else if(NR%4==2){print}}' ./reads/NC_000914_plasmid.fastq > ./reads/NC_000914_plasmid.fasta

awk '{if(NR%4==1){print ">chromosome"}else if(NR%4==2){print}}' ./reads/NC_000918_chromosome.fastq > ./reads/NC_000918_chromosome.fasta
awk '{if(NR%4==1){print ">chromosome"}else if(NR%4==2){print}}' ./reads/NC_002570_chromosome.fastq > ./reads/NC_002570_chromosome.fasta

cat ./reads/*.fasta > ./reads/reads.fasta
awk '{if(NR%2==1){print substr($1, 2)}}' ./reads/reads.fasta > ./reads/ids