/*
 * GEMMA 
 * Diego Garrido Martín 
 */

/*
 *  Define parameters
 */

// General params
params.genotype = 'data/mouse_hs1940.geno.txt'
params.phenotype = 'data/mouse_hs1940.pheno.txt'
params.pid = '1 2'
params.covariates = 'null'
params.dir = 'result'
params.out = 'gemma'
params.l = 7000
params.help = false

// Expand pid
if (params.pid =~ /,/) {
    pid = params.pid.replaceAll(',',' ')
} else if (params.pid =~ /:/) {
    def (start, end) = params.pid.tokenize(':')
    def seq = []
    def val = start.toInteger()
    while (val <= end.toInteger()) {
        seq += val
        val += 1
    }
    pid = seq.join(" ")
}


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
  log.info ' --genotype GENOTYPES        genotype data in GEMMA format (default: data/genotypes.vcf.gz)'
  log.info ' --phenotype PHENOTYPES      phenotype data in GEMMA format (default: data/phenotypes.tsv)'
  log.info ' --covariates COVARIATES     covariate data in GEMMA format (default: null)'
  log.info ' --l VARIANTS/CHUNK          number of variants tested per chunk (default: 10,000)'
  log.info ' --dir DIRECTORY             output directory (default: result)'
  log.info ' --out OUTPUT                output file prefix (default: gemma)'
  log.info ''
  exit(1)
}

/*
 *  Print parameter selection
 */

log.info ''
log.info 'Parameters'
log.info '------------------'
log.info "Genotype data                : ${params.genotype}"
log.info "Phenotype data               : ${params.phenotype}"
log.info "Covariates                   : ${params.covariates}"
log.info "Phenotype IDs                : ${params.pid}"
log.info "No. of variants/chunk        : ${params.l}"
log.info "Output directory             : ${params.dir}"
log.info "Output file prefix           : ${params.out}"
log.info ''


/*
 * Checks 
 */

// Mandatory options

if (!params.genotype) {
    exit 1, "Genotype file not specified."
} else if (!params.phenotype){
    exit 1, "Phenotype file not specified."
} 


/* 
 *  Split
 */

process split {

    input:
    file(geno) from file(params.genotype)

    output:
    file ("chunk*") into chunks_ch

    script:
    """
    split -d -a 10 -l ${params.l} $geno chunk
    """
}

/*
 *  Kinship
 */

process kinship {

    conda '/nfs/users2/rg/dgarrido/.conda/envs/gemma'

    input:
    file(geno) from file(params.genotype)
    file(pheno) from file(params.phenotype)

    output:
    file ('kinship.cXX.txt') into kinship_ch
    
    script:
    """
    gemma -g $geno -p $pheno -gk -outdir . -o kinship
    """
}


/*
 * Run gemma
 */

process gemma {

    conda '/nfs/users2/rg/dgarrido/.conda/envs/gemma'

    input:
    file(geno) from file(params.genotype)
    file(pheno) from file(params.phenotype)
    file(opt) from file(params.covariates)
    file(kinship) from kinship_ch 
    each file(chunk) from chunks_ch

    output:
    file ("${chunk}.assoc.out") into gemma_ch

    script:
    def cov = opt.name != 'null' ? "-c $opt" : '' 
    """
    gemma -g $geno -p $pheno -n $pid -k $kinship -lmm -outdir . -o $chunk $cov
    awk '{print \$2"\t"\$NF}' ${chunk}.assoc.txt > ${chunk}.assoc.out
    """
}


/* 
 * Reduce and generate output
 */

gemma_ch.collectFile(name: "${params.out}.assoc.txt", sort: { it.name }).set{out_ch}

process out {

    publishDir "${params.dir}"

    input:
    file(sim) from out_ch
 
    output:
    file(sim) into pub_ch
 
    script: 
    """
    
    """
}

