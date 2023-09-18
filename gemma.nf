/*
 * Copyright (c) 2021, Diego Garrido-Martín
 *
 * This file is part of 'gemma-nf':
 * A Nextflow pipeline for multivariate GWAS using GEMMA
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program. If not, see <http://www.gnu.org/licenses/>.
 */


// Define parameters

params.genotype = null
params.phenotype = null
params.maf = 0.01
params.dir = 'result'
params.out = 'gemma.tsv'
params.l = 500
params.t = 1
params.help = false


// Print usage

if (params.help) {
  log.info ''
  log.info 'G E M M A - N F'
  log.info '============================================================================================'
  log.info 'Performs multi-trait GWAS using GEMMA (https://github.com/genetics-statistics/GEMMA)'
  log.info ''
  log.info 'Usage: '
  log.info '    nextflow run gemma.nf [options]'
  log.info ''
  log.info 'Parameters:'
  log.info " --geno GENOTYPES            genotype data in VCF format, indexed"
  log.info " --pheno PHENOTYPES          (covariate-adjusted) phenotype data"
  log.info " --maf MAF                   MAF filter (default: ${params.maf}"
  log.info " --l VARIANTS/CHUNK          number of variants tested per chunk (default: ${params.l})"
  log.info " --t THREADS                 number of threads (deafault: ${params.t})"
  log.info " --dir DIRECTORY             output directory (default: ${params.dir})"
  log.info " --out OUTPUT                output file prefix (default: ${params.out})"
  log.info ''
  exit(1)
}


// Check mandatory parameters

if (!params.geno) {
    exit 1, "Genotype file not specified."
} else if (!params.pheno){
    exit 1, "Phenotype file not specified."
}


// Print parameter selection

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


// Split VCF

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


// Pre-process genotypes, phenotypes

process preprocess {

    cpus params.t

    input:
    file vcf from file(params.geno)    
    file pheno from file(params.pheno)

    output:
    tuple file("geno.bed"), file("geno.bim"), file("geno.fam") into geno_ch

    script:
    """
    comm -12 <(bcftools view -h $vcf | grep CHROM | sed 's/\\t/\\n/g' | sed '1,9d' | sort) <(cut -f1 $pheno | sort) > keep.txt
    plink2 --vcf $vcf --keep keep.txt --make-bed --out geno --threads ${params.t} 
    awk '{print int(NR)"\t"\$0}' <(cut -f2 geno.fam) > idx
    join -t \$'\t' -1 2 -2 1 <(sort -k2,2 idx) <(sort -k1,1 $pheno) | sort -k2,2n | cut -f1,2 --complement > pheno.tmp
    paste <(cut -f1-5 geno.fam) pheno.tmp > tmpfile; mv tmpfile geno.fam
    """
}


// Obtain kinship matrix

process kinship {

    cpus params.t

    input:
    tuple file(bed), file(bim), file(fam) from geno_ch

    output:
    tuple file("kinship.sXX.eigenD.txt"), ("kinship.sXX.eigenU.txt") into kinship_ch

    script:
    """
    # Compute kinship
    export OPENBLAS_NUM_THREADS=${params.t}
    prefix=\$(basename $bed | sed 's/.bed//')
    plink2 --bfile \$prefix --maf ${params.maf} --indep-pairwise 50 5 0.8 --threads ${params.t}
    plink2 --bfile \$prefix --extract plink2.prune.in --out geno.pruned --make-bed --threads ${params.t}
    gemma -gk 2 -bfile geno.pruned -outdir . -o kinship
    gemma -bfile geno.pruned -k kinship.sXX.txt -eigen -outdir . -o kinship.sXX
    """
}


// GWAS: testing (GEMMA)

process test {

    cpus params.t

    input:
    tuple file(bed),file(bim),file(fam) from geno_ch
    tuple file(kinship_d), file(kinship_u) from kinship_ch
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
            plink2 -bfile geno --extract bed1 region --make-bed --out geno.ss --threads ${params.t}
            paste <(cut -f1-5 geno.ss.fam) <(cut -f1-5 --complement geno.fam) > tmpfile; mv tmpfile geno.ss.fam
            (timeout 120 gemma -lmm -b geno.ss -d $kinship_d -u $kinship_u -n \$pids -outdir . -o gemma.k\$k -maf ${params.maf} &> STATUS || exit 0)
            if [[ \$(grep ERROR STATUS) ]]; then
                touch gemma.k\$k.assoc.txt
            else
                gemma -lmm -b geno.ss -d $kinship_d -u $kinship_u -n \$pids -outdir . -o gemma.k\$k -maf ${params.maf}
            fi
            ((k++))
        done
        cat gemma.k*.assoc.txt > gemma.\${chunknb}.assoc.txt        
    else
        paste <(head -1 $chunk) <(tail -1 $chunk | cut -f2) > region
        plink2 -bfile geno --extract bed1 region --make-bed --out geno.ss --threads ${params.t}
        paste <(cut -f1-5 geno.ss.fam) <(cut -f1-5 --complement geno.fam) > tmpfile; mv tmpfile geno.ss.fam
        (timeout 120 gemma -lmm -b geno.ss -d $kinship_d -u $kinship_u -n \$pids -outdir . -o gemma.\${chunknb} -maf ${params.maf} &> STATUS || exit 0)
        if [[ \$(grep ERROR STATUS) ]]; then
            touch gemma.\${chunknb}.assoc.txt
        else
            gemma -lmm -b geno.ss -d $kinship_d -u $kinship_u -n \$pids -outdir . -o gemma.\${chunknb} -maf ${params.maf}
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
