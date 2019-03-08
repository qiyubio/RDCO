for i in {1..19} X; do sh /home/qiyu/data/hiC/juicebox dump eigenvector KR ${1} chr${i} chr${i} bp 500000 ${1}_${i}.out; done

for i in {1..19} X; do cat ${1}_${i}.out >>${1}_merge.out; done

paste mm10_500k.txt_noY ${1}_merge.out|grep -v NaN > ${1}_merge.out.EI
