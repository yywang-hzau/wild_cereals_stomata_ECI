# Set SENTIEON_LICENSE if it is not set in the environment
module load SAMtools/1.9
module load sentieon/201808.07
export SENTIEON_LICENSE=mn01:9000

# Update with the location of the Sentieon software package

release_dir=/public/home/software/opt/bio/software/Sentieon/201808.07
#记得bwa要索引
# Other settings
nt=16 

#Sentieon proprietary compression
bam_option="--bam_compression 1"

# speed up memory allocation malloc in bwa
export LD_PRELOAD=$SENTIEON_INSTALL_DIR/lib/libjemalloc.so
export MALLOC_CONF=lg_dirty_mult:-1

#reference
fasta='td_genome.fa'
#i=$1
i=CRR0616692

fq1=${i}_1.paired.re.fq.gz
fq2=${i}_2.paired.re.fq.gz

group_prefix="ec1"
platform="ILLUMINA"
mq=30

# 输出文件
#rawCram=$i.cram

sortedCram=$i.q$mq.sorted.cram
sortedCram=CRR061692.sorted.cram
depCram=$i.deduped.cram
realnCram=$i.realn.cram
outvcf=$i.vcf
exec > $i.callVCF.log 2>&1

# ******************************************
# 1. BWA-MEM 
# ******************************************

#( $release_dir/bin/sentieon bwa mem -M -R "@RG\tID:${i}\tSM:${i}\tPL:$platform" \
#-t $nt -K 10000000 $fasta $fq1 $fq2 || echo -n 'error' ) | samtools sort -@ $nt  --output-fmt CRAM \
#--reference $fasta -o $rawCram - && samtools index -@ $nt $rawCram 
#samtools view -hCS -T $fasta -q $mq -o $sortedCram $rawCram && \
samtools index -@ $nt $sortedCram
samtools flagstat $rawCram > $i.stat.raw.txt && \
samtools flagstat $sortedCram > $i.stat.q$mq.txt &

# ******************************************
# 2. Calculate data metrics
# ******************************************
$release_dir/bin/sentieon driver -r $fasta -t $nt -i $sortedCram --algo MeanQualityByCycle ${i}_mq_metrics.txt \
--algo QualDistribution ${i}_qd_metrics.txt --algo GCBias --summary ${i}_gc_summary.txt ${i}_gc_metrics.txt \
--algo AlignmentStat --adapter_seq '' ${i}_aln_metrics.txt --algo InsertSizeMetricAlgo ${i}_is_metrics.txt
$release_dir/bin/sentieon plot metrics -o ${i}_metrics-report.pdf gc=${i}_gc_metrics.txt \
qd=${i}_qd_metrics.txt mq=${i}_mq_metrics.txt isize=${i}_is_metrics.txt
$release_dir/bin/sentieon driver -r $fasta -t $nt -i $sortedCram --algo LocusCollector --fun score_info ${i}_score.txt

# ******************************************
# 3. remove Duplicate Reads
# ******************************************
$release_dir/bin/sentieon driver -r $fasta -t $nt -i $sortedCram --algo Dedup --rmdup --cram_write_options version=3.0 \
--score_info ${i}_score.txt --metrics ${i}_dedup_metrics.txt $depCram 

# ******************************************
# 4. Indel 
# ******************************************
$release_dir/bin/sentieon driver -r $fasta -t $nt -i $depCram --algo Realigner  --cram_write_options version=3.0 \
$realnCram 

# ******************************************
# 5. Variant calling
# ******************************************
$release_dir/bin/sentieon driver -t $nt -r $fasta -i $realnCram --algo Genotyper  ${outvcf}
