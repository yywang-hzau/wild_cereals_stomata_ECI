vcftools --vcf  td_filter2.vcf.recode.vcf --chr chr1A --recode --recode-INFO-all --out td_chr1A.vcf
vcftools --vcf  td_filter2.vcf.recode.vcf --chr chr1B --recode --recode-INFO-all --out td_chr1B.vcf

vcftools --vcf  td_filter2.vcf.recode.vcf --chr chr2A --recode --recode-INFO-all --out td_chr2A.vcf
vcftools --vcf  td_filter2.vcf.recode.vcf --chr chr2B --recode --recode-INFO-all --out td_chr2B.vcf

vcftools --vcf  td_filter2.vcf.recode.vcf --chr chr3A --recode --recode-INFO-all --out td_chr3A.vcf
vcftools --vcf  td_filter2.vcf.recode.vcf --chr chr3B --recode --recode-INFO-all --out td_chr3B.vcf

vcftools --vcf  td_filter2.vcf.recode.vcf --chr chr4A --recode --recode-INFO-all --out td_chr4A.vcf
vcftools --vcf  td_filter2.vcf.recode.vcf --chr chr4B --recode --recode-INFO-all --out td_chr4B.vcf

vcftools --vcf  td_filter2.vcf.recode.vcf --chr chr5A --recode --recode-INFO-all --out td_chr5A.vcf
vcftools --vcf  td_filter2.vcf.recode.vcf --chr chr5B --recode --recode-INFO-all --out td_chr5B.vcf

vcftools --vcf  td_filter2.vcf.recode.vcf --chr chr6A --recode --recode-INFO-all --out td_chr6A.vcf
vcftools --vcf  td_filter2.vcf.recode.vcf --chr chr6B --recode --recode-INFO-all --out td_chr6B.vcf

vcftools --vcf  td_filter2.vcf.recode.vcf --chr chr7A --recode --recode-INFO-all --out td_chr7A.vcf
vcftools --vcf  td_filter2.vcf.recode.vcf --chr chr7B --recode --recode-INFO-all --out td_chr7B.vcf

cat td_chr1A.vcf.recode.vcf | awk '($1=="chr1A"){print $1":"$2,"chr1A", $2/100000000,$2,$4,$5}' >  chr1A.snp
cat td_chr2A.vcf.recode.vcf | awk '($1=="chr2A"){print $1":"$2,"chr2A", $2/100000000,$2,$4,$5}' >  chr2A.snp
cat td_chr3A.vcf.recode.vcf | awk '($1=="chr3A"){print $1":"$2,"chr3A", $2/100000000,$2,$4,$5}' >  chr3A.snp
cat td_chr4A.vcf.recode.vcf | awk '($1=="chr4A"){print $1":"$2,"chr4A", $2/100000000,$2,$4,$5}' >  chr4A.snp
cat td_chr5A.vcf.recode.vcf | awk '($1=="chr5A"){print $1":"$2,"chr5A", $2/100000000,$2,$4,$5}' >  chr5A.snp
cat td_chr6A.vcf.recode.vcf | awk '($1=="chr6A"){print $1":"$2,"chr6A", $2/100000000,$2,$4,$5}' >  chr6A.snp
cat td_chr7A.vcf.recode.vcf | awk '($1=="chr7A"){print $1":"$2,"chr7A", $2/100000000,$2,$4,$5}' >  chr7A.snp

cat td_chr1B.vcf.recode.vcf | awk '($1=="chr1B"){print $1":"$2,"chr1B", $2/100000000,$2,$4,$5}' >  chr1B.snp
cat td_chr2B.vcf.recode.vcf | awk '($1=="chr2B"){print $1":"$2,"chr2B", $2/100000000,$2,$4,$5}' >  chr2B.snp
cat td_chr3B.vcf.recode.vcf | awk '($1=="chr3B"){print $1":"$2,"chr3B", $2/100000000,$2,$4,$5}' >  chr3B.snp
cat td_chr4B.vcf.recode.vcf | awk '($1=="chr4B"){print $1":"$2,"chr4B", $2/100000000,$2,$4,$5}' >  chr4B.snp
cat td_chr5B.vcf.recode.vcf | awk '($1=="chr5B"){print $1":"$2,"chr5B", $2/100000000,$2,$4,$5}' >  chr5B.snp
cat td_chr6B.vcf.recode.vcf | awk '($1=="chr6B"){print $1":"$2,"chr6B", $2/100000000,$2,$4,$5}' >  chr6B.snp
cat td_chr7B.vcf.recode.vcf | awk '($1=="chr7B"){print $1":"$2,"chr7B", $2/100000000,$2,$4,$5}' >  chr7B.snp
