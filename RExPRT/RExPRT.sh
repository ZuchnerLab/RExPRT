#!/bin/bash
set -eou pipefail

if [ "$1" = "-h" -o "$1" = "--help" ]
then
    echo "    RExPRT is a machine learning tool to predict tandem repeat pathogenicity

    Usage: To run RExPRT use the following command:
        ./RExPRT.sh TRfile.txt
    TR file should have 5 columns labelled: chr, start, end, motif, and sampleID

    Example of a test file input is located in ./example/input/
    Example of test output files are in ./example/output/

    For full documentation and detailed instructions visit https://github.com/ZuchnerLab/RExPRT"
    exit 0
fi

if [ $# -eq 0 ]; then
    echo "No repeats file provided"
    exit 1
fi

#Download Annotation files
if [ ! -d "./annotation_files" ]
then
    wget -qO- https://zuchnerlab.s3.amazonaws.com/RExPRT_public/annotation_files.tar.gz | tar xvz
fi

#Download GERP files
if [ ! -d "./gerp_files" ]
then
    wget -qO- https://zuchnerlab.s3.amazonaws.com/RExPRT_public/gerp_files.tar.gz | tar xvz
fi

# Download helper scripts
if [ ! -d "./helper_scripts" ]
then
    wget -qO- https://zuchnerlab.s3.amazonaws.com/RExPRT_public/helper_scripts.tar.gz | tar xvz
fi


# Download S2SNet
if [ ! -d "./S2SNet" ]
then
    wget -qO- https://zuchnerlab.s3.amazonaws.com/RExPRT_public/S2SNet.tar.gz | tar xvz
fi

#Download ML models
if [ ! -f "SVM.pckl" ]
then
    wget https://zuchnerlab.s3.amazonaws.com/RExPRT_public/SVM.pckl
fi

if [ ! -f "XGB.pckl" ]
then
    wget https://zuchnerlab.s3.amazonaws.com/RExPRT_public/XGB.pckl
fi

# download bedtools
if [ ! -f "bedtools.static.binary" ]
then
    wget -q https://github.com/arq5x/bedtools2/releases/download/v2.29.2/bedtools.static.binary
    chmod u+x bedtools.static.binary
fi

# Annotate variants
repeats=$1

head -1 $1 > header
tail -n +2 $1 > repeats
./bedtools.static.binary sort -i repeats > del && mv del repeats
cat header repeats > sorted_repeats


echo "Annotating Exon and Intron"
head -q -n 1 sorted_repeats annotation_files/Exons_and_introns_UCSC.sorted.bed | paste -sd "\t" > header
awk '{print $0 "\tgene_distance"}' header > head && mv head header
sed 's/#//g' header > changed.txt && mv changed.txt header
./bedtools.static.binary closest -a sorted_repeats -b annotation_files/Exons_and_introns_UCSC.sorted.bed -d > intersection 
Rscript --vanilla helper_scripts/exon_intron_ann.R intersection header annotation_files/UCSC_canonical.txt
echo "Finished Exon Intron annotation"

chmod u+x helper_scripts/gerp_ann.sh
chmod u+x helper_scripts/split_bed.sh
echo "Annotating gerp scores"
awk '{print $1}' sorted_repeats | uniq | grep -wv "chr" > list.txt
time ./helper_scripts/gerp_ann.sh
echo " Finished annotating gerp scores"


echo "Annotating TAD boundaries"
head -1 final_annotated.txt > header
awk '{print $0 "\tTAD"}' header > head && mv head header
./bedtools.static.binary intersect -a final_annotated.txt -b annotation_files/TADboundaries_CpGcount.bed -c > intersection
cat header intersection > final_annotated.txt
echo "Finished annotating TAD boundaries"

echo "Annotating eSTR"
head -1 final_annotated.txt > header
awk '{print $0 "\teSTR"}' header > head && mv head header
./bedtools.static.binary intersect -a final_annotated.txt -b annotation_files/eSTR_loci_hg19.sorted.bed -c > intersection
cat header intersection > final_annotated.txt
echo "Finished annotating eSTR"

echo "Annotating opRegRegions"
head -1 final_annotated.txt > header
awk '{print $0 "\topReg"}' header > head && mv head header
./bedtools.static.binary intersect -a final_annotated.txt -b annotation_files/openRegulatoryRegions_hg19.sorted.bed -c > intersection
cat header intersection > final_annotated.txt
echo "Finished annotating opRegRegions"

echo "Annotating GTEx"
Rscript --vanilla helper_scripts/gtex_ann.R final_annotated.txt annotation_files/max_tissueExpression_perGene.txt
echo "Finished annotating GTEX"

echo "Annotating promoter regions"
head -1 final_annotated.txt > header
awk '{print $0 "\tpromoter"}' header > head && mv head header
./bedtools.static.binary intersect -a final_annotated.txt -b annotation_files/promoters.sorted.bed -c > intersection
cat header intersection > final_annotated.txt
echo "Finished annotating promoter regions"

echo "Annotating 3'UTR"
head -1 final_annotated.txt > header
awk '{print $0 "\tUTR_3"}' header > head && mv head header
./bedtools.static.binary intersect -a final_annotated.txt -b annotation_files/3primeUTR.sorted.bed -c > intersection
cat header intersection > final_annotated.txt
echo "Finished annotating 3'UTR"

echo "Annotating 5'UTR"
head -1 final_annotated.txt > header
awk '{print $0 "\tUTR_5"}' header > head && mv head header
./bedtools.static.binary intersect -a final_annotated.txt -b annotation_files/5primeUTR.sorted.bed -c > intersection
cat header intersection > final_annotated.txt
echo "Finished annotating 5'UTR"

echo "Annotating pLi and loeuf scores"
head -1 final_annotated.txt > header
awk '{print $0 "\tchrom\tcstart\tcend\tloeuf\tpLi"}' header > head && mv head header
./bedtools.static.binary intersect -a final_annotated.txt -b annotation_files/pLI_scores_hg19.sorted.bed -loj > intersection
cat header intersection > final_annotated.txt
Rscript --vanilla helper_scripts/pLi_ann.R final_annotated.txt
echo "Finished annotating pLi and loeuf scores"

echo "Annotating RAD21 binding sites"
head -1 final_annotated.txt > header
awk '{print $0 "\tRAD21"}' header > head && mv head header
cut -f2- annotation_files/NeuralCell_RAD21bindingSites_hg19.txt > file
./bedtools.static.binary intersect -a final_annotated.txt -b file -c > intersection
cat header intersection > final_annotated.txt
echo "Finished annotating RAD21"

echo "Annotating SMC3"
head -1 final_annotated.txt > header
awk '{print $0 "\tSMC3"}' header > head && mv head header
cut -f2- annotation_files/NeuralCell_SMC3bindingSites_hg19.txt > file
./bedtools.static.binary intersect -a final_annotated.txt -b file -c > intersection
cat header intersection > final_annotated.txt
echo "Finished annotating SMC3"

echo "Annotating percent GC"
Rscript --vanilla helper_scripts/perGC.R final_annotated.txt
echo "Finished Annotating percent GC"

echo "Adding S2S annotations"
pip install numpy --quiet
Rscript --vanilla helper_scripts/create_S2S_files.R final_annotated.txt
(cd S2SNet || exit 1; python S2SNet_noGUI_emb_py3_repeats.py) > text_delete
cut --complement -d$'\t' -f1,2,3 S2SNet/S2SNetTIs_Emb.txt > values
paste final_annotated.txt values > combined && mv combined final_annotated.txt
echo "Finished annotating S2S markers"


echo "Converting columns into binary"
Rscript --vanilla helper_scripts/convert_cols_binary.R final_annotated.txt
echo "Finished converting columns into binary"

echo "format for machine learning"
Rscript helper_scripts/format_final_annotated.R final_annotated.txt
echo "Done"


rm repeats
rm file
rm header
rm intersection
rm text_delete
rm values
rm list.txt
rm -r repeats_by_chrom
rm -r intersections
rm -r gerp_annotated/
rm combined_gerp.txt
rm sorted_repeats

#Use ML models to score annotated repeats
pip install sklearn -q
pip install xgboost -q
pip install category_encoders -q

python3 helper_scripts/rexprt.py

Rscript --vanilla helper_scripts/remove_duplicates.R TRsAnnotated_RExPRTscoresDups.txt RExPRT_scoresDups.txt

base=$(basename "$repeats")
filename="${base%.*}"
mv TRsAnnotated_RExPRTscores.txt ${filename}_TRsAnnotated_RExPRTscores.txt
mv RExPRTscores.txt ${filename}_RExPRTscores.txt
rm final_annotated.txt
rm TRsAnnotated_RExPRTscoresDups.txt
rm RExPRT_scoresDups.txt
