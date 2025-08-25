
/*
 * Pipeline for running the MolEpi GWAS scripts, microbiota branch optimized for large number of phenotypes.
 *
 * Per Lundmark, Molecular Epidemiology, Uppsala University
 * 
 * Based in part on shell code from Uwe Menzel, Molecular Epidemiology, Uppsala University
 *
 */

nextflow.enable.dsl=2

//Process chromosome list
clean_chr = params.chr.replaceAll("\\s","") //Remove any whitespace
chr_entries = clean_chr.split(',')
chromosomes = []
chr_entries.each{
    clean_chr_entry = it.replaceAll("\\s","")
    curr_entry = clean_chr_entry.split('-')
    if (curr_entry.size() == 1){
        chromosomes.addAll(curr_entry)
    }else if(curr_entry.size() == 2){
        chromosomes.addAll(curr_entry[0]..curr_entry[1])
    }else{
        println("Unable to parse chromosomes")
        system.exit(1)
    }
}
//Channel with specified chromosomes for analysis
chr_channel = Channel.fromList(chromosomes)


//Process phenotype name list
if (params.phenoname != "all"){
    clean_phenoname = params.phenoname.replaceAll("\\s","") //Remove any whitespace
    pn_entries = clean_phenoname.split(',')
    phenonames = []
    pn_entries.each{
        phenonames.add(it)
    }
    phenoname_channel = Channel.fromList(phenonames)
}

//Channel with phenotypes
pheno_channel = Channel.fromPath(params.phenofile)


/*
 * Log init, etc. before regression,
 * moved here from master gwas shell script (run_gwas.sh)
 */

process init_log{
    label 'short_job_1_core'
    publishDir "./${params.id}/association" //Collect output into this dir

    output:
    path '*_gwas.log'

    script:

    """
    log="${params.id}_gwas.log"   # master logfile 
    echo ""  > \${log}
    echo ""   | tee -a \${log}
    echo -n "  "  | tee -a \${log}
    date | tee -a \${log}
    echo "  Account: \${params.account}" | tee -a \${log}
    echo -n "  Operated by: " | tee -a \${log} 
    whoami | tee -a \${log} 
    echo "  Job identifier: " ${params.id} | tee -a \${log}
    echo "  Master logfile: \${log}" | tee -a \${log}
    echo "" | tee -a \${log}
    echo "  Running on chromosomes ${chromosomes}" | tee -a \${log}  
    echo "  Genotype input folder: ${params.genofolder}"  | tee -a \${log}
    echo "  Genotype identifier: ${params.genoid}" | tee -a \${log}
    echo "  Phenotype input file: ${params.phenofile}"  | tee -a \${log} 
    echo "  Phenotype column name(s): ${params.phenoname}"  | tee -a \${log}
    echo "  Covariate input file: ${params.covarfile}"  | tee -a \${log} 
    echo "  Covariate column name(s): ${params.covarnames}"  | tee -a \${log}
    echo "" | tee -a \${log}
    echo "  Threshold for minor allele count (mac): ${params.mac}" | tee -a \${log}
    echo "  Threshold for minor allele frequency (maf): ${params.maf}" | tee -a \${log}
    echo "  Maximum variance inflation factor (vif): ${params.vif}" | tee -a \${log}
    echo "  Maximum missing call rate for markers: ${params.geno} " | tee -a \${log}
    echo "  Maximum missing call rate for samples: ${params.mind} " | tee -a \${log}
    echo "  Threshold for Hardy-Weinberg p-value: ${params.hwe}"  | tee -a \${log}
    echo "  Imputation quality range (mach-r2): ${params.mr2}" | tee -a \${log}
    echo "" | tee -a \${log}
    echo "  Requested partition: ${params.gwas_partition}" | tee -a \${log}
    echo "" | tee -a \${log}
    """
}



process list_phenotypes{
    label 'short_job_1_core'
    
    output:
    path "phenotype_names.list"

    """
    head -n1 ${params.phenofile} | tr '\t| ' '\n' | tail -n +3 > phenotype_names.list
    """
}



process check_pheno_nonmissing{
    label 'short_job_1_core'
    label 'python3'
    publishDir "${params.id}/association"

    input:
    val phenotype_name

    output:
    tuple val(phenotype_name), path("nonmissing_*.list")

    script:
    if (params.pheno_specific_batch && params.covarnames){
        combined_covars = params.covarnames + ",batch_" + phenotype_name
    }
    if (params.covarnames){
        combined_covars = params.covarnames
    }
    if (params.pheno_specific_batch){
        combined_covars = "batch_" + phenotype_name
    }

    """
    python3 ${params.script_dir}/list_no_NA_samples.py -p ${params.phenofile} -pn ${phenotype_name} -c ${params.covarfile} -cn ${combined_covars} -o nonmissing_${phenotype_name}.list
    """
}



process run_regenie_step1{
    label 'long_job_16_core'
    label 'plink2'
    publishDir "${params.id}/association"

    input:
    tuple val(phenotype_name), path(nonmiss_list)

    output:
    tuple val(phenotype_name), path('fit_out*')  //, val(converged)

    script:
    step1_pgen_prefix = "${params.genofolder}/${params.regenie_step1_geno}"

    //Check if covars set for options to allow skipping if unset
    if (params.covarnames){
        covarnames_checked = "--covarColList " + params.covarnames
    }else{
        covarnames_checked = ""
    }

    
    if (params.covarnames_cat){
        cat_covars = params.covarnames_cat
        if (params.pheno_specific_batch){
            cat_covars = cat_covars + ",batch_" + phenotype_name
        }
        covarnames_cat_checked = "--catCovarList " + cat_covars
    }else if (params.pheno_specific_batch){
        covarnames_cat_checked = "--catCovarList batch_" + phenotype_name
    }else{
        covarnames_cat_checked = ""
    }
        

    """
    #Skipping logging for now, pendig removal or not from plink2 implementation, should slim down shell code and adapt groovy solutions

    #Filter on maf per phenotype to avoid low variance error in analysis
    plink2 \
        --pfile ${step1_pgen_prefix} \
        --keep ${nonmiss_list} \
        --maf 0.01 \
        --write-snplist \
        --out perpheno_maf0.01_${phenotype_name}

    regenie \
        --step 1 \
        --pgen ${step1_pgen_prefix} \
        --keep ${nonmiss_list} \
        --extract perpheno_maf0.01_${phenotype_name}.snplist \
        --covarFile ${params.covarfile} ${covarnames_checked} ${covarnames_cat_checked} \
        --phenoFile ${params.phenofile} --phenoColList ${phenotype_name} \
        --${params.trait_type} \
        --bsize ${params.bsize1} \
        --lowmem \
        --lowmem-prefix tmp_rg \
        --out fit_out_${phenotype_name} \
        --maxCatLevels 200 \
        --par-region b37 \
        --niter 300 \
        ${params.universal_param}

    """

}



process list_nonconverged_phenotypes{
    label 'short_job_1_core'
    publishDir "${params.id}/association"

    input: 
    tuple val(phenotype_name), path(step1_fit_out)

    output:
    path('non-converged_phenotypes.txt')

    script:
    """
    echo $phenotype_name >  non-converged_phenotypes.txt
    """
}



process run_regenie_step2{
    label 'long_job_16_core'
    publishDir "${params.id}/association"
    errorStrategy 'ignore'

    input:
    each chr
    tuple val(phenotype_name), path(step1_fit_out)

    output:
    //path '*_gwas_chr*.log', emit: assoc_log
    //path '*.regenie', emit: assoc_results
    tuple val(phenotype_name), path('*.regenie.gz'), path('*_gwas_chr*.log')

    script:
    pgen_prefix = "${params.genofolder}/${params.genoid}_${chr}"
    logchr      = "${params.id}_gwas_chrom${chr}.log"

    //Check if covars set for options to allow skipping if unset
    if (params.covarnames){
        covarnames_checked = "--covarColList " + params.covarnames
    }else{
        covarnames_checked = ""
    }

    if (params.covarnames_cat || params.pheno_specific_batch){
        cat_covars = params.covarnames_cat
        if(params.pheno_specific_batch){
            cat_covars = cat_covars + ",batch_" + phenotype_name
        }
        covarnames_cat_checked = "--catCovarList " + cat_covars 
    }else{
        covarnames_cat_checked = ""
    }

    """
    #Skipping logging for now, pendig removal or not from plink2 implementation, should slim down shell code and adapt groovy solutions
    
    #--minMAC ${params.mac} \
    #--minINFO ${params.info} \


    #outfile_prefix="${params.id}_gwas_chr"${chr} # OBS!! : Naming convention also used elsewhere!
    outfile_prefix="${params.id}_gwas_chr${chr}_${phenotype_name}" # OBS!! : Naming convention also used elsewhere!

    regenie \
        --step 2 \
        --pgen ${pgen_prefix} \
        --covarFile ${params.covarfile} ${covarnames_checked} ${covarnames_cat_checked} \
        --phenoFile ${params.phenofile} --phenoColList ${phenotype_name} \
        --${params.trait_type} \
        --bsize ${params.bsize2} \
        --pred fit_out_${phenotype_name}_pred.list \
        --out \${outfile_prefix} \
        --gz \
        --maxCatLevels 200 \
        --par-region b37 \
        --minMAC ${params.mac} \
        --minINFO ${params.info} \
        --niter 300 \
        ${params.universal_param}
    """
}



//Convert association test output to cojo input format
process convert_cojo_input_regenie{
    label 'std_job_1_core'
    label 'python3'
    publishDir "${params.id}/cojo"

    input:
    tuple val(phenotype_name), path(assoc_results), path(assoc_log)

    output:
    //tuple val(phenotype_name), path(assoc_results), path(assoc_log)
    tuple val(phenotype_name), path('*.ma')

    script:
    """
    #Combine regenie output for all chr into one file
    tail -n +2 -q *.regenie > sumstat_${phenotype_name}.regenie
    
    #Convert to cojo input format
    python3 ${params.script_dir}/convert_regenie_to_cojo.py -i sumstat_${phenotype_name}.regenie -o ${params.id}_cojo_input_${phenotype_name}.ma

    """

}

//Run cojo selection of independent associated variants
process run_cojo {
    label 'long_job_16_core'
    label 'gcta'

    input:
    each chr
    tuple val(phenotype_name), path(cojo_input_file)

    output:
    tuple val(phenotype_name), path('*.jma.cojo') 


    script:

    """
    touch ${params.id}_${phenotype_name}_chr${chr}.jma.cojo
    gcta64  --bfile ${params.genofolder}/${params.cojo_ref_bed_genoid}_${chr}  --thread-num 16 --chr ${chr} --cojo-file ${cojo_input_file} --cojo-slct --cojo-p ${params.cojo_pval} --cojo-wind ${params.cojo_window} --cojo-collinear ${params.cojo_colline} --out ${params.id}_${phenotype_name}_chr${chr}
    """

//ref for dev from settings
//cojo settings
//params.cojo_pval = '5e-8'
//params.cojo_window = '10000'
//params.cojo_colline = '0.9'
//params.cojo_maf = '' //Need this? Should be set at higher level in the script not specific to the cojo analysis?
}



/*
 * Create a parameter file with information for the 
 * review script.
 */

process create_parameterfile{
    label 'short_job_1_core'
    publishDir "${params.id}/association"

    output:
    path '*_gwas_params.txt'

    script:

    """
    paramfile="${params.id}_gwas_params.txt"    # Naming convention is also used in "review-GWAS.R". Do NOT change here without changing there!  
    echo -n > \${paramfile}
    workfolder="${params.id}"
    echo "plink2_version ${params.plink2_version}" >> \${paramfile}
    echo "workfolder \${workfolder}" >> \${paramfile} 
    echo "ident ${params.id}" >> \${paramfile}
    echo "cstart ${chromosomes[0]}" >> \${paramfile} 
    echo "cstop ${chromosomes[-1]}" >> \${paramfile}
    echo "genotype_id ${params.genoid}" >> \${paramfile}
    echo "phenofile ${params.phenofile}" >> \${paramfile}
    echo "phenoname ${params.phenoname}" >> \${paramfile}
    echo "covarfile ${params.covarfile}" >> \${paramfile}
    echo "covarname ${params.covarnames}" >> \${paramfile}
    echo "mac ${params.mac}" >> \${paramfile}
    echo "maf ${params.maf}" >> \${paramfile}
    echo "vif ${params.vif}" >> \${paramfile}
    echo "sample_max_miss ${params.mind}" >> \${paramfile}
    echo "marker_max_miss ${params.geno}" >> \${paramfile}
    echo "hwe_pval ${params.hwe}" >> \${paramfile} 
    machr2_low=\$( echo ${params.mr2}|awk '{print \$1}' )
    machr2_high=\$( echo ${params.mr2}|awk '{print \$2}' )
    echo "machr2_low \${machr2_low}" >> \${paramfile}
    echo "machr2_high \${machr2_high}" >> \${paramfile}

    echo "regenie_bsize1 ${params.bsize1}" >> \${paramfile}
    echo "regenie_bsize2 ${params.bsize2}" >> \${paramfile}
    echo "regenie_trait_type ${params.trait_type}" >> \${paramfile}
    echo "regenie_info ${params.info}" >> \${paramfile}
    echo "regenie_covarnames_cat ${params.covarnames_cat}" >> \${paramfile}
    """
}



/*
 * Finish up logging 
 */

process finish_log{
    label 'short_job_1_core'
    publishDir "${params.id}/association"

    input:
    path 'old.log' 

    output:
    path '*_gwas.log'

    script:
    newlog = params.id + "_gwas.log"
    //println(newlog)

    """
    cp old.log ${newlog}
    echo "" | tee -a ${newlog}
    echo -n "  "  | tee -a ${newlog}
    date | tee -a ${newlog}
    echo "  Done." | tee -a ${newlog}
    echo "" | tee -a ${newlog}
    """
}



process create_report_regenie {
    label 'long_job_16_core'
    label 'R_4_1_1'
    publishDir "${params.id}/report"

    input:
    path gwas_params
    //val pname
    val chromosomes
    //path assoc_results
    //path assoc_log

    tuple val(pname), path(assoc_results), path(assoc_log)
    //tuple val(phenotype_name), path('*.regenie'), path('*_gwas_chr*.log')

    output:
    path '*.html'
    path '*.RData'
    path 'lambda_*.txt', emit: lambda
    
    

    
    script:
    chrom_str = chromosomes.join(',')
    
    """
    #echo "${chrom_str}"

    ${params.script_dir}/review_gwas_regenie.R ${params.id} ${pname} ${chrom_str} ${params.script_dir}



    """
}



process concat_lambdas{
    label 'short_job_1_core'
    publishDir "${params.id}/report"

    input:
    path single_lambda 

    output:
    path "lambdas.tsv"

    script:
    """
    cat *.txt > lambdas.tsv
    """


}



// Main workflow definition
workflow{
    init_log()
    
    if (params.assoc == "regenie"){
        if (params.phenoname == "all"){
            list_phenotypes()
            check_pheno_nonmissing(list_phenotypes.out.splitText().map{it -> it.trim()})
        }else{
            check_pheno_nonmissing(phenoname_channel)
        }
        run_regenie_step1(check_pheno_nonmissing.out)
        //run_regenie_step1.out.view()

        step1_branched = run_regenie_step1.out.branch
        {
            converged: (it[1].size() == 3)
            non_converged: (it[1].size() == 2)
            error_state: true
        }
        //step1_branched.converged.view()
        //step1_branched.non_converged.view()

        run_regenie_step2(chr_channel, step1_branched.converged)
        list_nonconverged_phenotypes(step1_branched.non_converged.collect())

        if(params.run_cojo){
            convert_cojo_input_regenie(run_regenie_step2.out.groupTuple())
            run_cojo(chr_channel, convert_cojo_input_regenie.out.groupTuple())
        }
    }

    create_parameterfile()
    finish_log(init_log.out)

    if(params.assoc == 'regenie'){
        create_report_regenie(create_parameterfile.out, chr_channel.collect(), run_regenie_step2.out.groupTuple())
        concat_lambdas(create_report_regenie.out.lambda.collect())
    }
}


//Add check for consequences of variation for example the ensembl tool instead of just gene link and nearest gene stuff?
