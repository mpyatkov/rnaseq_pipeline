#! /bin/bash
########################################################################################
#Andy Rampersaud, 09.21.2015
########################################################################################
#This script is used to perform overlap between 2 BED files
#Sample command:
#./overlap.sh G80_M1_MACS_peaks.bed G80_M2_MACS_peaks.bed
if [ $# -ne 2 ]
then
  echo "Usage: `basename $0` peak1.bed peak2.bed"
  exit 
fi
# echo $1
# echo $2
peak1_base=`basename $1`
#echo $peak1_base
peak2_base=`basename $2`
#echo $peak2_base
peak1_name=${peak1_base%\.bed};
echo $peak1_name
peak2_name=${peak2_base%\.bed};
echo $peak2_name
########################################################################################
#BED3 doesn't work with intersectBed commands (need at least 9 columns for intersectBed commands to work)
#Otherwise I get the following error:
#Error: Type checker found wrong number of fields while tokenizing data line.
#Input BED files are reformatted to BED9 format
########################################################################################
echo "----------------------------------------------------------------------------"
echo "cleaning input"
awk 'BEGIN {OFS="\t"}
     {if ($1~/chr/ && $1 !="chrM" && $1 !~/random/ && $3>$2 && $2>0 && $3>0)
          print $1"\t"$2"\t"$3"\t""'$peak1_name'""\t"1000"\t"".""\t"0"\t"0"\t""0,0,0" > "peak1.bed";
      else 
          print $0 > "peak1_dump.bed"}' $1 
awk 'BEGIN {OFS="\t"}
     {if ($1~/chr/ && $1 !="chrM"  && $1 !~/random/ && $3>$2  && $2>0 && $3>0)
          print $1"\t"$2"\t"$3"\t""'$peak2_name'""\t"1000"\t"".""\t"0"\t"0"\t""0,0,0" > "peak2.bed";
      else 
          print $0 > "peak2_dump.bed"}' $2 
###############################
#Need a scratch dir since creating multiple temp files/removing files
Output=$peak1_name"_"$peak2_name"_Output"
#Need if statement when run script mulitple times
if [ ! -d $Output ]; then
mkdir $Output
else
rm $Output/*;
fi
cp peak1.bed ./$Output
cp peak2.bed ./$Output
rm peak1.bed
rm peak2.bed
cd $Output
###############################
echo "sorting input files"
sort -k1,1 -k2,2n peak1.bed > peak1.temp
mv peak1.temp peak1.bed
sort -k1,1 -k2,2n peak2.bed > peak2.temp
mv peak2.temp peak2.bed
###############################
echo "classify common or unique peaks"
intersectBed -a peak1.bed -b peak2.bed -u | sort -k1,1 -k2,2n -k3,3n > $peak1_name'_common.bed'
intersectBed -a peak2.bed -b peak1.bed -u | sort -k1,1 -k2,2n -k3,3n > $peak2_name'_common.bed'
intersectBed -a peak1.bed -b peak2.bed -v | sort -k1,1 -k2,2n -k3,3n > $peak1_name'_unique.bed'
intersectBed -a peak2.bed -b peak1.bed -v | sort -k1,1 -k2,2n -k3,3n > $peak2_name'_unique.bed'
echo "getting the union of the peak sets"
cat peak1.bed peak2.bed > temp1.bed
sort -k1,1 -k2,2n temp1.bed > temp2.bed
bedtools merge -i temp2.bed > temp3.bed
cp temp3.bed $peak1_name'_'$peak2_name'_Union.bed'
mv temp3.bed Union_Peak_Set.bed
rm temp*.bed
echo "getting the merged common peaks"
cat $peak1_name'_common.bed' $peak2_name'_common.bed' > temp1.bed
sort -k1,1 -k2,2n temp1.bed > temp2.bed
bedtools merge -i temp2.bed > $peak1_name'_'$peak2_name'_merged_common.bed'
cp $peak1_name'_'$peak2_name'_merged_common.bed' Merged_Common_Peak_Set.bed
rm temp*.bed
echo "getting distance of nearest peak"
closestBed -d -a peak1.bed -b peak2.bed > $peak1_name'_distance.bed'
closestBed -d -a peak2.bed -b peak1.bed > $peak2_name'_distance.bed'
awk '{print $NF}' $peak1_name'_distance.bed' > temp1.txt
mv temp1.txt $peak1_name'_distance.txt'
rm $peak1_name'_distance.bed'
awk '{print $NF}' $peak2_name'_distance.bed' > temp2.txt
mv temp2.txt $peak2_name'_distance.txt'
rm $peak2_name'_distance.bed'
echo "counting peaks"
peak1_total=$(wc -l < 'peak1.bed') 
peak2_total=$(wc -l < 'peak2.bed') 
peak1_common=$(wc -l < $peak1_name'_common.bed') 
peak2_common=$(wc -l < $peak2_name'_common.bed') 
peak1_unique=$(wc -l < $peak1_name'_unique.bed') 
peak2_unique=$(wc -l < $peak2_name'_unique.bed')
echo "Running overlap.R script to get peak overlap, peak width, and peak proximity statistics ..."
cd ..
Rscript overlap.R $peak1_total $peak2_total $peak1_common $peak2_common $peak1_unique $peak2_unique $peak1_name $peak2_name Union_Peak_Set Merged_Common_Peak_Set
#echo "Rscript overlap.R $peak1_total $peak2_total $peak1_common $peak2_common $peak1_unique $peak2_unique $peak1_name $peak2_name Union_Peak_Set Merged_Common_Peak_Set"
cd $Output
rm peak1*
rm peak2*
rm Union_Peak_Set.bed
rm Merged_Common_Peak_Set.bed
echo "Creating summary file"
cat *_Stats.txt >> Peak_Width.txt
cat *_Stats_2.txt >> Peak_Width_2.txt
cat *_distance_Summary.txt >> Distance_Summary.txt
###################################
#Add label to Distance_Summary.txt
echo 'Peak proximity(bp):' >> Header3.txt
cat Header3.txt Distance_Summary.txt >> temp3.txt
mv temp3.txt Distance_Summary.txt
rm Header3.txt
###################################
#Add label to Peak_Width.txt
echo 'Peak width(bp):' >> Header1.txt
cat Header1.txt Peak_Width.txt >> temp1.txt
mv temp1.txt Peak_Width.txt
rm Header1.txt
###################################
#Add label to Peak_Width_2.txt
echo 'Peak width(bp):' >> Header1.txt
cat Header1.txt Peak_Width_2.txt >> temp1.txt
mv temp1.txt Peak_Width_2.txt
rm Header1.txt
###################################
#Add label to *_Overlap_Summary.txt
echo 'Peak overlap:' >> Header2.txt
cat Header2.txt *_Overlap_Summary.txt >> temp2.txt
mv temp2.txt  *_Overlap_Summary.txt
rm Header2.txt
###################################
echo >> *_Overlap_Summary.txt
cat *_Overlap_Summary.txt Peak_Width.txt >> Overlap_Summary.temp
echo >> Overlap_Summary.temp
cat Overlap_Summary.temp Distance_Summary.txt >> Overlap_Summary.temp2
echo >> Overlap_Summary.temp2 
cat Overlap_Summary.temp2 Peak_Width_2.txt >> Overlap_Summary.txt
#Remove extra files
rm Overlap_Summary.temp
rm Overlap_Summary.temp2 
rm *_Stats.txt
rm *_Overlap_Summary.txt
rm Peak_Width.txt
rm *_distance.txt
rm *_distance_Summary.txt
rm Distance_Summary.txt
rm *_2.txt
cd ..
echo "Done!"
echo "----------------------------------------------------------------------------"
echo
