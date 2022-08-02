dir=//

for i in `ls *.tped`
do
i=${i/.tped/}
cut -d " " -f 5- ${i}.tped  | awk '{print $0" "}' > ${i}.geno
done

xpclr --format txt --out ${dir}chr1A_xp --map ${dir}chr1A.snp  --popA ${dir}Chr1A.popas.geno --popB ${dir}Chr1A.popes.geno --chr chr1A --phased  --size 2000000 --step 500000
xpclr --format txt --out ${dir}chr2A_xp --map ${dir}chr2A.snp  --popA ${dir}Chr2A.popas.geno --popB ${dir}Chr2A.popes.geno --chr chr2A --phased  --size 2000000 --step 500000
xpclr --format txt --out ${dir}chr3A_xp --map ${dir}chr3A.snp  --popA ${dir}Chr3A.popas.geno --popB ${dir}Chr3A.popes.geno --chr chr3A --phased  --size 2000000 --step 500000
xpclr --format txt --out ${dir}chr4A_xp --map ${dir}chr4A.snp  --popA ${dir}Chr4A.popas.geno --popB ${dir}Chr4A.popes.geno --chr chr4A --phased   --size 2000000 --step 500000
xpclr --format txt --out ${dir}chr5A_xp --map ${dir}chr5A.snp  --popA ${dir}Chr5A.popas.geno --popB ${dir}Chr5A.popes.geno --chr chr5A --phased  --size 2000000 --step 500000
xpclr --format txt --out ${dir}chr6A_xp --map ${dir}chr6A.snp  --popA ${dir}Chr6A.popas.geno --popB ${dir}Chr6A.popes.geno --chr chr6A --phased --size 2000000 --step 500000
xpclr --format txt --out ${dir}chr7A_xp --map ${dir}chr7A.snp  --popA ${dir}Chr7A.popas.geno --popB ${dir}Chr7A.popes.geno --chr chr7A --phased --size 2000000 --step 500000

xpclr --format txt --out ${dir}chr1B_xp --map ${dir}chr1B.snp  --popA ${dir}Chr1B.popas.geno --popB ${dir}Chr1B.popes.geno --chr chr1B --phased  --size 2000000 --step 500000
xpclr --format txt --out ${dir}chr2B_xp --map ${dir}chr2B.snp  --popA ${dir}Chr2B.popas.geno --popB ${dir}Chr2B.popes.geno --chr chr2B --phased  --size 2000000 --step 500000
xpclr --format txt --out ${dir}chr3B_xp --map ${dir}chr3B.snp  --popA ${dir}Chr3B.popas.geno --popB ${dir}Chr3B.popes.geno --chr chr3B --phased  --size 2000000 --step 500000
xpclr --format txt --out ${dir}chr4B_xp --map ${dir}chr4B.snp  --popA ${dir}Chr4B.popas.geno --popB ${dir}Chr4B.popes.geno --chr chr4B --phased   --size 2000000 --step 500000
xpclr --format txt --out ${dir}chr5B_xp --map ${dir}chr5B.snp  --popA ${dir}Chr5B.popas.geno --popB ${dir}Chr5B.popes.geno --chr chr5B --phased  --size 2000000 --step 500000
xpclr --format txt --out ${dir}chr6B_xp --map ${dir}chr6B.snp  --popA ${dir}Chr6B.popas.geno --popB ${dir}Chr6B.popes.geno --chr chr6B --phased --size 2000000 --step 500000
xpclr --format txt --out ${dir}chr7B_xp --map ${dir}chr7B.snp  --popA ${dir}Chr7B.popas.geno --popB ${dir}Chr7B.popes.geno --chr chr7B --phased --size 2000000 --step 500000
