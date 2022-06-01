#!/bin/bash                                                                                                                                          
# align reads to the genome

cd ../nohmrRNA/

array=(038_SG_Gale_HIV_Lat_IFNy_JLat92_8h_RNA_mRNA_run03_S2  039_SG_Gale_HIV_Lat_IFNy_JLat92_8h_RNA_mRNA_run03_S3    040_SG_Gale_HIV_Lat_IFNy_JLat92_12h_RNA_mRNA_run03_S4   041_SG_Gale_HIV_Lat_IFNy_JLat92_12h_RNA_mRNA_run03_S5   042_SG_Gale_HIV_Lat_IFNy_JLat92_12h_RNA_mRNA_run03_S6)


#array=(001_SG_Gale_HIV_Lat_Jurkat_4h_RNA_mRNA_merged	002_SG_Gale_HIV_Lat_Jurkat_4h_RNA_mRNA_merged	003_SG_Gale_HIV_Lat_Jurkat_4h_RNA_mRNA_run01_S3	004_SG_Gale_HIV_Lat_IFNB_Jurkat_4h_RNA_mRNA_merged	005_SG_Gale_HIV_Lat_IFNB_Jurkat_4h_RNA_mRNA_merged	006_SG_Gale_HIV_Lat_IFNB_Jurkat_4h_RNA_mRNA_merged	007_SG_Gale_HIV_Lat_IFNB_Jurkat_8h_RNA_mRNA_merged	008_SG_Gale_HIV_Lat_IFNB_Jurkat_8h_RNA_mRNA_merged	009_SG_Gale_HIV_Lat_IFNB_Jurkat_8h_RNA_mRNA_merged	010_SG_Gale_HIV_Lat_IFNB_Jurkat_12h_RNA_mRNA_merged	011_SG_Gale_HIV_Lat_IFNB_Jurkat_12h_RNA_mRNA_merged	012_SG_Gale_HIV_Lat_IFNB_Jurkat_12h_RNA_mRNA_merged	013_SG_Gale_HIV_Lat_IFNy_Jurkat_4h_RNA_mRNA_run01_S13	014_SG_Gale_HIV_Lat_IFNy_Jurkat_4h_RNA_mRNA_merged	015_SG_Gale_HIV_Lat_IFNy_Jurkat_4h_RNA_mRNA_run01_S15	016_SG_Gale_HIV_Lat_IFNy_Jurkat_8h_RNA_mRNA_run01_S16	017_SG_Gale_HIV_Lat_IFNy_Jurkat_8h_RNA_mRNA_run01_S17	018_SG_Gale_HIV_Lat_IFNy_Jurkat_8h_RNA_mRNA_merged	019_SG_Gale_HIV_Lat_IFNy_Jurkat_12h_RNA_mRNA_merged	020_SG_Gale_HIV_Lat_IFNy_Jurkat_12h_RNA_mRNA_run02_S2	021_SG_Gale_HIV_Lat_IFNy_Jurkat_12h_RNA_mRNA_run02_S3	022_SG_Gale_HIV_Lat_JLat92_4h_RNA_mRNA_merged	023_SG_Gale_HIV_Lat_JLat92_4h_RNA_mRNA_merged	024_SG_Gale_HIV_Lat_JLat92_4h_RNA_mRNA_merged	025_SG_Gale_HIV_Lat_IFNB_JLat92_4h_RNA_mRNA_merged	026_SG_Gale_HIV_Lat_IFNB_JLat92_4h_RNA_mRNA_merged	027_SG_Gale_HIV_Lat_IFNB_JLat92_4h_RNA_mRNA_run02_S9	028_SG_Gale_HIV_Lat_IFNB_JLat92_8h_RNA_mRNA_merged	029_SG_Gale_HIV_Lat_IFNB_JLat92_8h_RNA_mRNA_run02_S11	030_SG_Gale_HIV_Lat_IFNB_JLat92_8h_RNA_mRNA_run02_S12	031_SG_Gale_HIV_Lat_IFNB_JLat92_12h_RNA_mRNA_run02_S13	032_SG_Gale_HIV_Lat_IFNB_JLat92_12h_RNA_mRNA_merged	033_SG_Gale_HIV_Lat_IFNB_JLat92_12h_RNA_mRNA_run02_S15	034_SG_Gale_HIV_Lat_IFNy_JLat92_4h_RNA_mRNA_merged	035_SG_Gale_HIV_Lat_IFNy_JLat92_4h_RNA_mRNA_merged	036_SG_Gale_HIV_Lat_IFNy_JLat92_4h_RNA_mRNA_merged	037_SG_Gale_HIV_Lat_IFNy_JLat92_8h_RNA_mRNA_run03_S1	038_SG_Gale_HIV_Lat_IFNy_JLat92_8h_RNA_mRNA_run03_S2	039_SG_Gale_HIV_Lat_IFNy_JLat92_8h_RNA_mRNA_run03_S3	040_SG_Gale_HIV_Lat_IFNy_JLat92_12h_RNA_mRNA_run03_S4	041_SG_Gale_HIV_Lat_IFNy_JLat92_12h_RNA_mRNA_run03_S5	042_SG_Gale_HIV_Lat_IFNy_JLat92_12h_RNA_mRNA_run03_S6)


echo "Array size: ${#array[*]}"
echo "Array items:"

for item in ${array[*]}
do

    printf "   %s\n" $item

    # STAR
    cmd='/vol01/ngs_tools/mapper/STAR-STAR_2.4.0h1/source/STAR --genomeDir /vol01/genome/human/STAR --clip5pNbases 1 --clip3pNbases 24 --readFilesCommand zcat --outSAMtype BAM SortedByCoordinate --readFilesIn '$item'_nohmrRNA.fastq.1.gz '$item'_nohmrRNA.fastq.2.gz --outFileNamePrefix ../mapping/'$item' --runThreadN 50 1>>../logs/mapping/'$item'_mapping.log 2>&1'

    #echo $cmd
    #echo $cmd > ../logs/mapping/$item\_mapping.log
    #eval $cmd
    
    # move a copy of the STAR generated log to the log folder
    #cp ../mapping/$item\Log.final.out ../logs/mapping/


    cmd2='/vol01/ngs_tools/python/install/Python-2.7.3/bin/htseq-count --stranded=reverse --format=bam --mode=intersection-nonempty --idattr=gene_id ../mapping/'$item'Aligned.sortedByCoord.out.bam.sub.bam /vol01/genome/human/igenome/Homo_sapiens/NCBI/build37.1/Annotation/Genes/genes.gtf > ../counts/'$item'_counts.txt &'

    # move a copy of the STAR generated log to the log folder
    #cp ../mapping/$item\Log.final.out ../logs/mapping/


 #   echo $cmd2
    echo $cmd2 > ../logs/mapping/$item\_mapping.log
   eval $cmd2

    # move a copy of the STAR generated log to the log folder
    #cp ../mapping/$item\Log.final.out ../logs/mapping/


done
