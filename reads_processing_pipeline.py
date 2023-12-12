####### Working Directory and Project Specifics ############
WRKDIR = '/mnt/dfc_data2/project/chenwentao/project/switch/data_analysis/kraken_clincal_samples/'
readWD = '/mnt/dfc_data2/project/chenwentao/project/switch/data_analysis/kraken_clincal_samples/'
SAMPLES, = glob_wildcards(readWD + 'symlinks/{samp}_1.clean.fq.gz')

################   REFERENCE    ##############
REF = '/mnt/dfc_data2/project/chenwentao/project/switch/primers_design/targetgenome/Chlamydia_trachomatis_DUW-3CX_chromosome_complete_genome.fasta'
Krakendb = '/mnt/dfc_data2/project/chenwentao/dataset/kraken2/maxikraken2_1903_140GB'

##########################################################################################
rule all:
    input:
#        expand("1.kraken2_classfication/{samp}_output.txt", samp=SAMPLES),
#        expand("2.Chlamydia_only_reads/{samp}_chlamydiae_only_1.fq", samp=SAMPLES),
#        expand("3.trimmed_reads/{samp}_R1.PAIREDtrimmomatictrimmed.fq", samp=SAMPLES),
#        expand("4.downsampled_reads/{samp}_R1.downsampled.fq", samp=SAMPLES),
#        expand("5.alignment/{samp}.dedup.bam", samp=SAMPLES),
#        expand("6.variants/masked_snps_{samp}.vcf.gz", samp=SAMPLES),
        expand("7.pseudosequences/{samp}.pseudosequence.uncovMasked.fasta", samp=SAMPLES)

#####################################################################


rule kraken2_classfication:
    input:
        r1="symlinks/{samp}_1.clean.fq.gz",
        r2="symlinks/{samp}_2.clean.fq.gz"
    output:
        out_txt="1.kraken2_classfication/{samp}_output.txt",
        report_txt="1.kraken2_classfication/{samp}_report.txt",
        classified_r1="1.kraken2_classfication/kraken_{samp}_1.fq",
        classified_r2="1.kraken2_classfication/kraken_{samp}_2.fq"
    shell:
        """
        kraken2 --db {Krakendb} --threads 12 --paired --classified-out 1.kraken2_classfication/kraken_{wildcards.samp}#.fq {input.r1} {input.r2} --output {output.out_txt} --report {output.report_txt}
        """

rule krakentools_extract_chlamydiae_reads:
    input:
        r1="1.kraken2_classfication/kraken_{samp}_1.fq",
        r2="1.kraken2_classfication/kraken_{samp}_2.fq",
        out_txt="1.kraken2_classfication/{samp}_output.txt",
        report_txt="1.kraken2_classfication/{samp}_report.txt"
    output:
        out_r1="2.Chlamydia_only_reads/{samp}_chlamydia_only_1.fq",
        out_r2="2.Chlamydia_only_reads/{samp}_chlamydia_only_2.fq"
    shell:
        """
        extract_kraken_reads.py -k {input.out_txt} -s1 {input.r1} -s2 {input.r2} -t 810 --include-children -o {output.out_r1} -o2 {output.out_r2} --report {input.report_txt} --fastq-output
        """

rule trim_illumina_Adaptors_fastqs:
    input:
        r1="2.Chlamydia_only_reads/{samp}_chlamydia_only_1.fq",
        r2="2.Chlamydia_only_reads/{samp}_chlamydia_only_2.fq"
    output:
        out_r1_paired="3.trimmed_reads/{samp}_R1.PAIREDtrimmomatictrimmed.fq",
        out_r1_unpaired="3.trimmed_reads/{samp}_R1.UNPAIREDtrimmomatictrimmed.fq",
        out_r2_paired="3.trimmed_reads/{samp}_R2.PAIREDtrimmomatictrimmed.fq",
        out_r2_unpaired="3.trimmed_reads/{samp}_R2.UNPAIREDtrimmomatictrimmed.fq"
    shell:
        """
        # Use -phred64 for simulated reads, -phred33 for NovaSeq 
        trimmomatic PE -phred64 {input.r1} {input.r2} {output.out_r1_paired} {output.out_r1_unpaired} {output.out_r2_paired} {output.out_r2_unpaired} ILLUMINACLIP:/mnt/dfc_data1/home/chenwentao/anaconda2/envs/kraken2/share/trimmomatic/adapters/TruSeq3-PE.fa:2:30:10:2:True LEADING:3 TRAILING:3 MINLEN:36
        """ 

rule seqtk_downsample:
    input:
        r1="3.trimmed_reads/{samp}_R1.PAIREDtrimmomatictrimmed.fq",
        r2="3.trimmed_reads/{samp}_R2.PAIREDtrimmomatictrimmed.fq"
    output:
        r1_downsampled="4.downsampled_reads/{samp}_R1.downsampled.fq",
        r2_downsampled="4.downsampled_reads/{samp}_R2.downsampled.fq"
    shell:
        """
        seqtk sample -s 100 {input.r1} 5000000 > {output.r1_downsampled}
        seqtk sample -s 100 {input.r2} 5000000 > {output.r2_downsampled}
        """

rule bwa_mem_mapping:
    input:
        ref=REF,
        r1="4.downsampled_reads/{samp}_R1.downsampled.fq",
        r2="4.downsampled_reads/{samp}_R2.downsampled.fq"
    output:
        bam="5.alignment/{samp}.bam"
    shell:
        """
        bwa mem -t 12 {input.ref} {input.r1} {input.r2} -Y -M -R "@RG\\tID:bwa\\tPL:illumina\\tLB:{wildcards.samp}_lib\\tSM:{wildcards.samp}" |
        samtools view -Sb - > {output.bam}
        """

rule samtools_sort_index:
    input:
        "5.alignment/{samp}.bam"
    output:
        bam="5.alignment/{samp}.sorted.bam",
        bai="5.alignment/{samp}.sorted.bam.bai"
    shell:
        """
        samtools sort {input} -o {output.bam}
        samtools index {output.bam}
        """

rule gatk_indel_realigner:
    input:
        bam="5.alignment/{samp}.sorted.bam",
        ref=REF
    output:
        intervals="5.alignment/{samp}.intervals",
        realigned_bam="5.alignment/{samp}.realigned.bam"
    shell:
        """
        gatk3 -T RealignerTargetCreator -R {input.ref} -I {input.bam} -o {output.intervals}
        gatk3 -T IndelRealigner -R {input.ref} -I {input.bam} -targetIntervals {output.intervals} -o {output.realigned_bam}
        """

rule picard_markduplicates:
    input:
        "5.alignment/{samp}.realigned.bam"
    output:
        dedup_bam="5.alignment/{samp}.dedup.bam",
        metrics="5.alignment/{samp}.metrics.txt"
    shell:
        """
        picard MarkDuplicates I={input} O={output.dedup_bam} M={output.metrics} REMOVE_DUPLICATES=true
        """
rule variant_calling:
    input:
        bam="5.alignment/{samp}.dedup.bam",
        ref=REF
    output:
        vcf="6.variants/{samp}.vcf.gz",
    shell:
        """
        bcftools mpileup -f {input.ref} -Ou {input.bam} | bcftools call --ploidy 1 -mv -Oz -o {output.vcf}
        """

rule extract_snps:
    input:
        vcf="6.variants/{samp}.vcf.gz"
    output:
        snps_vcf="6.variants/{samp}.snps.vcf.gz"
    shell:
        """
        bcftools view -v snps {input.vcf} -Oz -o {output.snps_vcf}
        bcftools index {output.snps_vcf}
        """

rule mask_variants:
    input:
        snps_vcf="6.variants/{samp}.snps.vcf.gz"
    output:
        masked_vcf="6.variants/masked_snps_{samp}.vcf.gz"
    shell:
        """
        bcftools view -i 'MQ<20 || DP<5 || INFO/DP4[2]<2 || INFO/DP4[3]<2 || ((INFO/DP4[2] + INFO/DP4[3]) / (INFO/DP4[0] + INFO/DP4[1] + INFO/DP4[2] + INFO/DP4[3])) < 0.8' {input.snps_vcf} -Oz -o {output.masked_vcf}
        """

rule generate_pseudosequence:
    input:
        snps_vcf="6.variants/{samp}.snps.vcf.gz",
        masked_vcf="6.variants/masked_snps_{samp}.vcf.gz",
        ref=REF
    output:
        pseudo_seq="7.pseudosequences/{samp}.filtered.pseudosequence.fasta"
    shell:
        """
        bcftools consensus -f {input.ref} {input.snps_vcf} -m {input.masked_vcf} -o {output.pseudo_seq}.tmp --mask-with N && \
        echo ">{wildcards.samp}" > {output.pseudo_seq} && \
        tail -n +2 {output.pseudo_seq}.tmp >> {output.pseudo_seq} && \
        rm {output.pseudo_seq}.tmp
        """

rule mask_uncovered_regions_in_pseudosequence:
    input:
        bam="5.alignment/{samp}.dedup.bam",
        pseudo_seq_filtered="7.pseudosequences/{samp}.filtered.pseudosequence.fasta"
    output:
        pseudo_seq="7.pseudosequences/{samp}.pseudosequence.uncovMasked.fasta"
    shell:
        """
        # Generate a BED file of uncovered regions
        samtools depth -a {input.bam} | awk '$3 == 0' | \
        awk -v samp="{wildcards.samp}" 'BEGIN{{OFS="\t"}} {{print samp, $2-1, $2}}' > 7.pseudosequences/{wildcards.samp}.uncovered.bed

        # Mask the uncovered regions in the pseudo-sequence using bedtools
        bedtools maskfasta -fi {input.pseudo_seq_filtered} -bed 7.pseudosequences/{wildcards.samp}.uncovered.bed -fo {output.pseudo_seq}
        """

