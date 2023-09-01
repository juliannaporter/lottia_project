# delly + manta files were .gz for bcftools merge (to create pop files)
# had to be unzipped for survivor (at least I think they do)
rule gunzip:
    input:
        delly = "data/delly/delly.split.wholegenome.sample/{sample}.vcf.gz",
        manta = "data/manta/manta.split.wholegenome.sample/{sample}.vcf.gz"
    output:
        delly = "data/delly/delly.split.wholegenome.sample/{sample}.vcf",
        manta = "data/manta/manta.split.wholegenome.sample/{sample}.vcf"
    shell:
        """
        module load htslib
        bgzip -d {input.delly}
        bgzip -d {input.manta}
        """

rule surv_callset:
    input:
        delly = ancient("data/delly/delly.split.wholegenome.sample/{sample}.vcf"), 
        manta = ancient("data/manta/manta.split.wholegenome.sample/{sample}.vcf")
    output:
        "data/survivor/{sample}.survivor.callset"
    shell:
        """
        echo {input.delly} > {output}
        echo {input.manta} >> {output}
        """

rule survivor:
    input:
        ancient("data/survivor/{sample}.survivor.callset")
    output:
        "data/survivor/{sample}.survivor.vcf"
    shell:
        """
        install/SURVIVOR/Debug/SURVIVOR merge \
        {input} 1000 2 1 0 0 50 {output}
        """ 

rule match_gt:
    input:
        "data/survivor/{sample}.survivor.vcf"
    output:
        delly = "data/GT_filt/{sample}.col10.dellyGT.txt",
        manta = "data/GT_filt/{sample}.col11.mantaGT.txt",
        GTmatch = "data/GT_filt/{sample}.matching.GT.txt",
        sitelist = "data/GT_filt/{sample}.GTmatch.sitelist.bed"
    shell:
        """
        module load bcftools
        bcftools view -H {input} | awk '{{print $1,$2, substr($11,0,3);}}' > {output.manta}
        bcftools view -H {input} | awk '{{print $1,$2, substr($10,0,3);}}' > {output.delly}
        grep -Fxf {output.manta} {output.delly} > {output.GTmatch}
        cut -d " " -f 1,2 {output.GTmatch} > {output.sitelist}
        """

rule filt_prep:
    input:
        vcf = "data/survivor/{sample}.survivor.vcf",
        region = "data/GT_filt/{sample}.GTmatch.sitelist.bed"
    output:
        vcf = "data/survivor/{sample}.survivor.sorted.vcf.gz",
        region = "data/GT_filt/{sample}.final.sitelist"
    params:
        vcf = "data/survivor/{sample}.survivor.sorted.vcf"
    shell:
        """
        sed 's# #\t#g' {input.region} > {output.region}
        module load bcftools
        bcftools sort -o {params.vcf} \
        {input.vcf}
        module load htslib
        bgzip {params.vcf}
        tabix {output.vcf}
        """

## ycuky stuff about this rule:
#### I had to remake the .bed file b/c the sep was a space not a tab
#### I had to sort, bgzip, and tabix the files
#### all was put into rule above
rule filt_gt:
    input:
        vcf = "data/survivor/{sample}.survivor.sorted.vcf.gz",
        region = "data/GT_filt/{sample}.final.sitelist"
    output:
        filt = "data/GT_filt/{sample}.gtfilt.vcf",
        col = "data/GT_filt/{sample}.gtfilt.1col.vcf"
    params:
        col = "{sample}"
    shell:
        """
        module load bcftools
        bcftools view \
        --regions-file {input.region} \
        -Ov -o {output.filt} \
        {input.vcf}
        bcftools view \
        --samples {params.col} \
        -Ov -o {output.col} \
        {output.filt}
        """

rule merge_prep:
    input:
        "data/GT_filt/{sample}.gtfilt.1col.vcf"
    output:
        "data/GT_filt/{sample}.gtfilt.1col.vcf.gz"
    shell:
        """
        module load htslib
        bgzip {input}
        tabix {output}
        """

rule check_ref:
    input:
        vcf = "data/GT_filt/{sample}.gtfilt.1col.vcf.gz",
        ref = config["ref"]
    output:
        vcf = "data/GT_filt/{sample}.gtfilt.1col.norm.vcf.gz",
        tbi = "data/GT_filt/{sample}.gtfilt.1col.norm.vcf.gz.tbi"
    shell:
        """
        module load bcftools
        bcftools norm --check-ref ws -Oz \
        -o {output.vcf} --fasta-ref {input.ref} {input.vcf}
        module load htslib
        tabix {output.vcf}
        """

rule merge_geno:
    input:
        expand("data/GT_filt/{sample}.gtfilt.1col.norm.vcf.gz.tbi", sample = SAMPLE)
    output:
        "data/final_cnv/All.samples.gtfilt.vcf.gz"
    shell:
        """
        module load bcftools
        bcftools merge \
        -Oz -o {output} \
        `ls data/GT_filt/*gtfilt.1col.norm.vcf.gz`
        """

### filtering final CNV callset 

rule max_no_call_loose:
    input:
        "data/final_cnv/All.samples.gtfilt.vcf.gz"
    output:
        "data/final_cnv/All.samples.gtfilt.maxnocall50.vcf.gz"
    shell:
        """
        module load GATK
        gatk SelectVariants \
        -V {input} \
        --exclude-filtered true \
        --max-nocall-fraction 0.50 \
        -O {output}
        """

rule max_no_call_strict:
    input:
        "data/final_cnv/All.samples.gtfilt.vcf.gz"
    output:
        "data/final_cnv/All.samples.gtfilt.maxnocall20.vcf.gz"
    shell:
        """
        module load GATK
        gatk SelectVariants \
        -V {input} \
        --exclude-filtered true \
        --max-nocall-fraction 0.20 \
        -O {output}
        """

rule maf:
    input:
        "data/final_cnv/All.samples.gtfilt.maxnocall{nocall}.vcf.gz"
    output:
        "data/final_cnv/All.samples.gtfilt.maxnocall{nocall}.maf.recode.vcf.gz"
    params:
        prefix = "data/final_cnv/All.samples.gtfilt.maxnocall{nocall}.maf"
    shell:
        """
        module load vcftools
        vcftools \
        --gzvcf {input} \
        --maf 0.05 \
        --remove-filtered-all \
        --recode \
        --out {params.prefix}
        """
