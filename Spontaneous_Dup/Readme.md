
### Code to detect spontaneous duplicates using Delly (https://github.com/dellytools/delly, Version: 0.8.1)


<br />


## Pipelines


#### Step 1) delly call on each bam

        find /data/neurospora/bam_files/ -name "*.bam"  -print | sort | \
            xargs -n 1 -P 4 -I BAM_FILE \
            sh -c '
                sample=`basename BAM_FILE | cut -d"." -f1`
                
                delly call \
                    -g reference_genome/neurospora_crassa.fasta \
                    -o /data/delly/${sample}.delly.bcf \
                    BAM_FILE
            '


#### Step 2) merge all results :

        BCF_FILES=`find /data/delly/ -name "*.delly.bcf" | \
            sort | xargs -I BCF echo -n "BCF " | sed 's/ $//'`
        delly merge -o /data/delly/all_samples.delly.bcf \
            ${BCF_FILES}


#### Step 3) re-genotyping on all detected positions

        find /data/neurospora/bam_files/ -name "*.bam"  -print | sort | \
            xargs -n 1 -P 4 -I BAM_FILE \
            sh -c '
                sample=`basename BAM_FILE | cut -d"." -f1`
                
                delly call \
                    -g reference_genome/neurospora_crassa.fasta \
                    -v /data/delly/all_samples.delly.bcf \
                    -o /data/delly/${sample}.delly.rec.bcf \
                    BAM_FILE
            '


#### Step 4) merge final delly results

        BCF_FILES=`find /data/delly/ -name "*.delly.rec.bcf" | \
            sort | xargs -I BCF echo -n "BCF " | sed 's/ $//'`
            
        bcftools merge -m id -O b \
            -o /data/delly/all_samples.delly.rec.bcf \
            ${BCF_FILES}


#### Step 5) filtering delly calls, and detect spontaneous duplicates which only present in asci against parents with a 2:2 segregation 

        ## Criteria used to filter Delly results (reference to https://www.biostars.org/p/116628/):
        ## PASS: >= 3 reads support with mapping quality (QUAL >=20)
        ## PRECISE: support by split reads


        ## Using Cross E as an example
        bcftools view --samples-file /data/Cross_E.samples.txt -c 1 \
            /data/delly/all_samples.delly.rec.bcf | \
            detect_mutations.pl -v - --max-cmp-depth 2 --max-cmp-total 3 --min-supp-depth 3 -a \
            -g /data/Cross_E.2_groups.txt | \
            awk '/\#/ || (!/GRPID=parent/ && /SVTYPE=DUP/ && $7 !~ /LowQual/)' \
            > /data/Cross_E.delly.candidate_dups.vcf


####   

    

<br />


