hic=$1
res=$2
chr=$3
title=$4
python /home/qiyu/data/hiC/HiCPlotter-master/HiCPlotter.py -f "$hic"_"$chr"_"$res".out_matrix.txt -chr "$chr" -r "$res" -o "$hic"_"$chr"_"$res" -ptr 1 -hist  "$hic"_"$chr"_"$res".out_matrix_withhead.is500001.ids200001.insulation.bedGraph.bedGraph,"$hic"_"$chr"_"$res".out_matrix_withhead.is500001.ids1000001.insulation.bedGraph.bedGraph,"$hic"_"$chr"_"$res".out_matrix_withhead.is500001.ids600001.insulation.bedGraph.bedGraph -hl r500k_w200k,r500k_w600k,r500k_w1m  -t "$hic"_"$chr"_"$res".out_matrix_withhead.is500001.ids200001.insulation.boundaries.bed.bed -tl r500k_w200k_boundary -n "$title"

#sh plot_figure.sh inter_30_preLep.hic 100000 chr15 PreLep
#sh plot_figure.sh inter_30_Lep.hic 100000 chr15 Lep
#sh plot_figure.sh inter_30_4N_ZPD.hic 100000 chr15 4N_ZPD
#sh plot_figure.sh WhT_B6XCAST_MboI_inter_30.hic 100000 chr15 Whole_Testis

