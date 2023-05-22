/*
 * GEMMA 
 * Diego Garrido Martín 
 */

/*
 *  Define parameters
 */

// General params
params.genotype = null
params.phenotype = null
params.maf = 0.01
params.dir = 'result'
params.out = 'gemma.tsv'
params.l = 100000
params.t = 1
params.help = false

/*
 *  Print usage and help
 */

if (params.help) {
  log.info ''
  log.info 'G E M M A - N F'
  log.info '======================================================================='
  log.info 'Genome-wide Efficient Mixed-Model Analysis for Association Studies'
  log.info ''
  log.info 'Usage: '
  log.info '    nextflow run gemma.nf [options]'
  log.info ''
  log.info 'Parameters:'
  log.info " --geno GENOTYPES            genotype data in PLINK format (default: ${params.geno})"
  log.info " --pheno PHENOTYPES          phenotype data (default: ${params.pheno})"
  log.info " --maf MAF                   MAF filter (default: ${params.maf}"
  log.info " --l VARIANTS/CHUNK          number of variants tested per chunk (default: ${params.l})"
  log.info " --t THREADS                 number of threads (deafault: ${params.t})"
  log.info " --dir DIRECTORY             output directory (default: ${params.dir})"
  log.info " --out OUTPUT                output file prefix (default: ${params.out})"
  log.info ''
  exit(1)
}

/*
 *  Print parameter selection
 */

log.info ''
log.info 'Parameters'
log.info '------------------'
log.info "Genotype data                : ${params.geno}"
log.info "Phenotype data               : ${params.pheno}"
log.info "MAF                          : ${params.maf}"
log.info "No. of variants/chunk        : ${params.l}"
log.info "No. of threads               : ${params.t}"
log.info "Output directory             : ${params.dir}"
log.info "Output file prefix           : ${params.out}"
log.info ''


/*
 * Checks 
 */

// Mandatory options

if (!params.geno) {
    exit 1, "Genotype file not specified."
} else if (!params.pheno){
    exit 1, "Phenotype file not specified."
}


/*
 *  Split
 */

process split {

    input:

    file(vcf) from file(params.geno)
    file(index) from file("${params.geno}.tbi")

    output:

    file("chunk*") into chunks_ch

    script:
    """
    bcftools query -f '%CHROM\t%POS\n' $vcf > positions
    split -d -a 10 -l ${params.l} positions chunk
    """
}

/*
 *  Prepare
 */

process preprocess {

    input:
    file vcf from file(params.geno)    
    file pheno from file(params.pheno)

    output:
    tuple file("geno.bed"), file("geno.bim"), file("geno.fam") into geno_ch

    script:
    """
    comm -12 <(bcftools view -h $vcf | grep CHROM | sed 's/\\t/\\n/g' | sed '1,9d' | sort) <(cut -f1 $pheno | sort) > keep.txt
    plink2 --vcf $vcf --keep keep.txt --make-bed --out geno
    awk '{print int(NR)"\t"\$0}' <(cut -f2 geno.fam) > idx
    join -t \$'\t' -1 2 -2 1 <(sort -k2,2 idx) <(sort $pheno) | sort -k2,2n | cut -f1,2 --complement > pheno.tmp
    paste <(cut -f1-5 geno.fam) pheno.tmp > tmpfile; mv tmpfile geno.fam
    """
}

/*
 *  Kinship
 */

process kinship {

    cpus params.t

    input:
    tuple file(bed), file(bim), file(fam) from geno_ch

    output:
    file("kinship.sXX.txt.gz") into kinship_ch

    script:
    """
    # Compute kinship
    export OPENBLAS_NUM_THREADS=${params.t}
    gemma -gk 2 -bfile \$(basename $bed | sed 's/.bed//') -outdir . -o kinship
    gzip kinship.sXX.txt
    """
}

/*
 *  Test
 */

process test {

    cpus params.t

    input:
    tuple file(bed),file(bim),file(fam) from geno_ch
    file(kinship) from kinship_ch
    each file(chunk) from chunks_ch

    output:
    file('gemma.0*.assoc.txt') into sstats_ch

    script:
    """
    export OPENBLAS_NUM_THREADS=${params.t}
    pids=\$(seq \$(echo \$(awk '{print NF}' $fam | head -1) - 5 | bc -l))
    chunknb=\$(basename $chunk | sed 's/chunk//')

    if [[ \$(cut -f1 $chunk | sort | uniq -c | wc -l) -ge 2 ]]; then
        k=1
        cut -f1 $chunk | sort | uniq | while read chr; do
            paste <(grep -P "^\$chr\t" $chunk | head -1) <(grep -P "^\$chr\t" $chunk | tail -1 | cut -f2) > region
            plink2 -bfile geno --extract bed1 region --make-bed --out geno.ss
            paste geno.ss.fam <(cut -f1-5 --complement geno.fam) > tmpfile; mv tmpfile geno.ss.fam
            (timeout 10 gemma -lmm -b geno.ss -k $kinship -n \$pids -outdir . -o gemma.k\$k -maf ${params.maf} &> STATUS || exit 0)
            if [[ \$(grep ERROR STATUS) ]]; then
                touch gemma.k\$k.assoc.txt
            else
                gemma -lmm -b geno.ss -k $kinship -n \$pids -outdir . -o gemma.k\$k -maf ${params.maf}
            fi
            ((k++))
        done
        cat gemma.k*.assoc.txt > gemma.\${chunknb}.assoc.txt        
    else
        paste <(head -1 $chunk) <(tail -1 $chunk | cut -f2) > region
        plink2 -bfile geno --extract bed1 region --make-bed --out geno.ss
        paste geno.ss.fam <(cut -f1-5 --complement geno.fam) > tmpfile; mv tmpfile geno.ss.fam
        (timeout 10 gemma -lmm -b geno.ss -k $kinship -n \$pids -outdir . -o gemma.\${chunknb} -maf ${params.maf} &> STATUS || exit 0)
        if [[ \$(grep ERROR STATUS) ]]; then
            touch gemma.\${chunknb}.assoc.txt
        else
            gemma -lmm -b geno.ss -k $kinship -n \$pids -outdir . -o gemma.\${chunknb} -maf ${params.maf}
        fi
    fi
    """
}


sstats_ch.collectFile(name: "${params.out}", sort: { it.name }).set{pub_ch}


// Summary stats

process end {

    publishDir "${params.dir}", mode: 'copy'

    input:
    file(out) from pub_ch

    output:
    file(out) into end_ch

    script:
    """
    head -1 $out > header
    grep -v n_miss $out > tmpfile; cat header tmpfile > $out
    """
}
