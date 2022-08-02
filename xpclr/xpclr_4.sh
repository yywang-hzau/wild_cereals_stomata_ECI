cat xpclr.txt | awk '{print$1,$2,$3,$12}' | grep -v 'e-' | sed 's/-[0-9.]/0/g' | awk  'IF$4>0{print$0}' | grep -v ' 0.0' | awk '{print$1"\t"$2"\t"$3"\t"$4}'| grep -v "chrom" | sort -t $'\t' -k 1,1 -k 2n,2 > xpclr_plot.txt
 cat xpclr_plot.txt | sort -k 4 -nr | head -n 21 | sort -t $'\t' -k 1,1 -k 2n,2> top5.txt
