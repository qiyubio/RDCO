cat $1|grep -v "e+"|perl -lane 'if ($F[1]>0 and $F[1]<$F[2]) {print $F[0],"\t",$F[1],"\t",$F[2],"\t",$F[3],"\t",0,"\t",$F[4]}' >${1}.bed
#grep "+" ${1}.bed >${1}.forward.bed
#grep "-" ${1}.bed >${1}.reverse.bed

#bedtools getfastq -fi ../genome/canFam3/genome.fa -bed ${1}.forward.bed -fo ${1}.forward.fa
#bedtools getfastq -fi ../genome/canFam3/genome.fa -bed ${1}.reverse.bed -fo ${1}.reverse.fa

#perl fasta_to_fastq.pl ${1}.forward.fa >${1}.forward.fq
#perl fasta_to_fastq.pl ${1}.reverse.fa >${1}.reverse.fq
#bedtools getfasta -fi ../genome/canFam3_flip/genome.fa -bed ${1}.bed -fo ${1}.fa
#perl fasta_to_fastq.pl ${1}.fa >${1}.fq
#bwa mem -t 24 /data/yuq3/genome/canFam3_flip/version0.7.10/genome.fa ${1}.fq|samtools view -bhS -|samtools sort - -o ${1}.sort.bam 
#bamToBed -i ${1}.sort.bam > ${1}.sortall.bed
#grep -P ${2}"\t" ${1}.sortall.bed >${1}.sort.bed
#bedToBam -i ${1}.bed -g canFam3_new.genome1 >${1}.bam
#samtools sort ${1}.bam >${1}.sort.bam
#samtools index ${1}.sort.bam
sort -k1,1 -k2,2n ${1}.bed >${1}.sort.bed
#grep -P ${2}"\t" canFam3_100kb_new.genome >canFam3_100kb_new.genome_${2}
coverageBed -b ${1}.sort.bed -a /data/yuq3/simulation/canFam3_100kb_new.genome >${1}_100kb.bedgraph
macs2 callpeak -t ${1}.sort.bed --bw=1000 --gsize=2100000000 -q=0.01 --slocal=2000 --llocal=10000 -m 3 30 --name=${1}_macs2 2>${1}_macs2.out

if test ! -e dog1_peak; then
cat ${2}|cut -f1-3|grep -v from|grep -v chrX >dog1_peak

fi


if test -f ${1}_macs2_peaks.xls; then 

#sh coverage_biowulf.sh ${1}.bed dog1_peak_${2} 4000

grep -vP '(^\#|chr\s|^$)' ${1}_macs2_peaks.xls|cut -f1,2,3,10 >${1}_macs2_peaks.xls.bed

intersectBed -a ${1}_macs2_peaks.xls.bed -b dog1_peak -u >${1}_overlap.bed
intersectBed -a ${1}_macs2_peaks.xls.bed -b dog1_peak -v >${1}_nonoverlap.bed
intersectBed -b ${1}_macs2_peaks.xls.bed -a dog1_peak -u >${1}_peak_at.bed

fi
