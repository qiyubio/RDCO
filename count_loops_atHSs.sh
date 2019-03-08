#!/bin/bash
## Author: Qi Yu
## Contact: dryuqi@gmail.com


cat ${1}|perl -lane '{print $F[0],"\t",$F[1],"\t",$F[2],"\t","hichip","\t",$F[7],"\t","+","\t",$F[4],"\t",$F[5],"\t","255,0,0","\t",2,"\t","20,20","\t","0,",$F[5]-$F[1]}' >${1}_convert



#cat "track name=${3}" ${1}_IGV.bed 

intersectBed -a ${1}_convert -b $2 -u >${1}_at_${2}

cat ${1}_at_${2}|perl -lane '{print $F[0],"\t",$F[6],"\t",$F[7],"\t",$F[3],"\t",$F[4],"\t",$F[5],"\t",$F[1],"\t",$F[2],"\t",$F[8],"\t",$F[9],"\t",$F[10],"\t",$F[11]}' >${1}_at_${2}_rev

intersectBed -a ${1}_at_${2}_rev -b $2 -u >${1}_at_${2}_both



cat ${1}_at_${2}_both|perl -lane ' if ($F[6]!=$F[1]|$F[7]!=$F[2]) {print $F[0],"\t",$F[6],"\t",$F[7],"\t",$F[3],"\t",$F[4],"\t",$F[5],"\t",$F[1],"\t",$F[2],"\t",$F[8],"\t",$F[9],"\t",$F[10],"\t",$F[11]}' >${1}_at_${2}_both_noself

cat ${1}_convert|perl -lane ' if ($F[6]!=$F[1]|$F[7]!=$F[2]) {print $F[0],"\t",$F[6],"\t",$F[7],"\t",$F[3],"\t",$F[4],"\t",$F[5],"\t",$F[1],"\t",$F[2],"\t",$F[8],"\t",$F[9],"\t",$F[10],"\t",$F[11]}' >${1}_noself


cat ${1}_convert|perl -lane '{print $F[0],"\t",$F[6],"\t",$F[7],"\t",$F[3],"\t",$F[4],"\t",$F[5],"\t",$F[1],"\t",$F[2],"\t",$F[8],"\t",$F[9],"\t",$F[10],"\t",$F[11]}' >${1}_convert_rev

intersectBed -a ${1}_convert_rev -b $2 -u|perl -lane 'if ($F[6]!=$F[1]|$F[7]!=$F[2]) {print $F[0],"\t",$F[6],"\t",$F[7],"\t",$F[3],"\t",$F[4],"\t",$F[5],"\t",$F[1],"\t",$F[2],"\t",$F[8],"\t",$F[9],"\t",$F[10],"\t",$F[11]}' >${1}_rev_at_${2}

cat ${1}_at_${2} ${1}_rev_at_${2}|sort|uniq|perl -lane 'if ($F[6]!=$F[1]|$F[7]!=$F[2]) {print $_}'|intersectBed -a - -b ${1}_at_${2}_both_noself -v >${1}_merge_at_${2}


cat ${1}_at_${2}_both_noself|perl -lane '{print $F[0],"\t",$F[1],"\t",$F[7],"\t",$F[3],"\t",$F[4],"\t",$F[5],"\t",$F[1],"\t",$F[7],"\t",$F[8],"\t",$F[9],"\t",$F[10],"\t",$F[11]}' >${1}_at_${2}_both_noself.IGV.bed

perl -lane 'if ($F[4]>=5) {print $_}' ${1}_at_${2}_both_noself.IGV.bed >${1}_at_${2}_both_noself_5.IGV.bed

perl -lane 'if ($F[4]>=5) {print $_}' ${1}_at_${2}_both_noself >${1}_at_${2}_both_noself5

overlap_hts_one5=$(wc -l <${1}_at_${2}_both_noself5);echo ${overlap_hts_one5}

sed -i '1s/^/track name="junctions"\n/' ${1}_at_${2}_both_noself.IGV.bed
sed -i '1s/^/track name="junctions"\n/' ${1}_at_${2}_both_noself_5.IGV.bed


overlap_hts=$(wc -l <${1}_at_${2}_both);echo ${overlap_hts}
total_loop=$(wc -l <${1}_convert);echo ${total_loop}

total_loop_nos=$(wc -l <${1}_noself);echo ${total_loop_nos}
overlap_hts_nos=$(wc -l <${1}_at_${2}_both_noself);echo ${overlap_hts_nos}

overlap_hts_one=$(wc -l <${1}_merge_at_${2});echo ${overlap_hts_one}


printf $1 >>loop_summary.txt

awk '''BEGIN {printf ("\t%s\t%s\t%.2f\t%s\t%s\t%.2f\t%s\t%.2f\t%s\t%.2f\n",'$total_loop','$overlap_hts','$overlap_hts'*100/'$total_loop','$total_loop_nos','$overlap_hts_nos','$overlap_hts_nos'*100/'$total_loop_nos',
'$overlap_hts_one5','$overlap_hts_one5'*100/'$total_loop_nos','$overlap_hts_one','$overlap_hts_one'*100/'$total_loop_nos')}''' >>loop_summary.txt





