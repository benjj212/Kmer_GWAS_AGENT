### comparing end of chromsome 2B between the different genomes
## first step is to extract the region of interest for each of genomes using bedtools
bedtools getfasta -fi genome.fasta -bed region_of_interest.bed > region_og_interest.fasta

## replacing the chromosome names in each of the files 
sed -i 's/>Taes_Norin61_v15_chr2B:782944975-802944975/>chr2B/g' Taes_Norin61_chr2B.fasta
sed -i 's/>Taes_Renan_v1_chr2B:792232696-812232696/>chr2B/g' Taes_Renan_chr2B.fasta
sed -i 's/>Chr2B:792755788-812755788/>chr2B/g' Chr2B_CS.fasta
sed -i 's/>Taes_Landmark_v1_v1_chr2B:777311483-797311483/>chr2B/g' Taes_Landmark_chr2B.fasta
sed -i 's/>Taes_Mace_v1_chr2B:776420183-796420183/>chr2B/g' Taes_Mace_chr2B.fasta
sed -i 's/>Taes_Stanley_v1_2_chr2B:770745243-790745243/>chr2B/g' Taes_Stanley_chr2B.fasta
sed -i 's/>Taes_SYMattis_v1_1_chr2B:779857935-799857935/>chr2B/g' Taes_SYMattis_chr2B.fasta
sed -i 's/>Taes_Lancer_chr2b.fasta:779857935-799857935/>chr2B/g' Taes_Lancer_chr2b.fasta

## aligning the genomes using Minimap
minimap2 -ax asm5 --eqx Chr2B_CS.fasta Taes_Landmark_chr2B.fasta > CS_LAND.sam
minimap2 -ax asm5 --eqx Taes_Landmark_chr2B.fasta Taes_Norin61_chr2B.fasta > LAND_NO.sam
minimap2 -ax asm5 --eqx Taes_Norin61_chr2B.fasta Taes_Stanley_chr2B.fasta > NO_STAN.sam
minimap2 -ax asm5 --eqx Taes_Stanley_chr2B.fasta Taes_Renan_chr2B.fasta > STAN_REN.sam
minimap2 -ax asm5 --eqx Taes_Renan_chr2B.fasta Taes_Mace_chr2B.fasta > REN_MACE.sam
minimap2 -ax asm5 --eqx Taes_Mace_chr2B.fasta Taes_SYMattis_chr2B.fasta > MACE_SYM.sam
minimap2 -ax asm5 --eqx Taes_SYMattis_chr2B.fasta Taes_Lancer_chr2b.fasta > SYM_LANC.sam

## extracting the alignment info using Syri
syri -c CS_LAND.sam -r Chr2B_CS.fasta -q Taes_Landmark_chr2B.fasta -k -F S --prefix CS_LAND
syri -c LAND_NO.sam -r Taes_Landmark_chr2B.fasta -q Taes_Norin61_chr2B.fasta -k -F S --prefix LAND_NO
syri -c NO_STAN.sam -r Taes_Norin61_chr2B.fasta -q Taes_Stanley_chr2B.fasta -k -F S --prefix NO_STAN
syri -c STAN_REN.sam -r Taes_Stanley_chr2B.fasta -q Taes_Renan_chr2B.fasta -k -F S --prefix STAN_REN
syri -c REN_MACE.sam -r Taes_Renan_chr2B.fasta -q Taes_Mace_chr2B.fasta -k -F S --prefix REN_MACE
syri -c MACE_SYM.sam -r Taes_Mace_chr2B.fasta -q Taes_SYMattis_chr2B.fasta -k -F S --prefix MACE_SYM
syri -c SYM_LANC.sam -r Taes_SYMattis_chr2B.fasta -q Taes_Lancer_chr2b.fasta -k -F S --prefix SYM_LANC

## plotting using the plotsr function 
plotsr --sr CS_LANDsyri.out \
       --sr LAND_NOsyri.out \
       --sr NO_STANsyri.out \
       --sr STAN_RENsyri.out \
       --sr REN_MACEsyri.out \
       --sr MACE_SYMsyri.out \
       --sr SYM_LANCsyri.out \
       --genomes genomes.txt -b pdf -H 20 -W 10 --markers markers.bed -o out.pdf



