### Generate synthetic (in silico) mutations from real sequencing reads. 
Those synthetic mutations could be used to test the false-negative rate of a certain mutation-detection pipeline.

<br />

#### Original methods described in:

* 1) Keightley, P. D., Ness, R. W., Halligan, D. L. & Haddrill, P. R. Estimation of the Spontaneous Mutation Rate per Nucleotide Site in a Drosophila melanogaster Full-Sib Family. Genetics 196, 313¨C320 (2014).

* 2) Keightley, P. D. et al. Estimation of the Spontaneous Mutation Rate in Heliconius melpomene. Mol Biol Evol 32, 239¨C243 (2015).

<br />



## Pipelines
* Using Cross E as an example, related files could be found under the "example_data/" folder

<br />

#### Step 1: get the empirical distributions of the depth of non-reference alleles

* This information could be simply extract from a VCF file or readcount tools. Note for hapoid genome, the non-reference allele ratio is 1 for nearly all sites except some with sequencing errors/mapping artefacts.

<br />

#### Step 2: generate synthetic reads (reads with synthetic mutations) for each tetrad/ascus

* Random select a genome position with a read depth x, then replace y reads with a random picked nucleotide (different from the original one), the y was determined according to the empirical distribution. Generate reads with synthesized mutations, excluding non-informative reads (-F 3844):    
	1. unmapped;   
	2. not primary alignment;   
	3. read fails platform/vendor quality checks;   
	4. read is PCR or optical duplicate;   
	5. supplementary alignment;      

<br />

* Generate 5,000 mutation sites with 2:2 ratio

        sim_mutation_reads.pl --fasta reference_genome/neurospora_crassa.fasta \
            --depth empirical_non-ref_allele_depths.tbl \
            --random-size 5000 --samtools "-F 3844" --exclude Supercontig_12.8 Supercontig_12.9 Supercontig_12.10 Supercontig_12.11 \
                Supercontig_12.12 Supercontig_12.13 Supercontig_12.14 Supercontig_12.15 Supercontig_12.16 Supercontig_12.17 \
                Supercontig_12.18 Supercontig_12.19 Supercontig_12.20 \
            --group-file Cross_E/Cross_E_s40.bam_groups.txt \
            --min-group-size 2 --max-group-size 2 \
            > Cross_E/Cross_E_s40.simulated.dat

* Generate 1,000 mutation sites with 3:1 ratio

        sim_mutation_reads.pl --fasta reference_genome/neurospora_crassa.fasta \
            --depth empirical_non-ref_allele_depths.tbl \
            --random-size 1000 --samtools "-F 3844" --exclude Supercontig_12.8 Supercontig_12.9 Supercontig_12.10 Supercontig_12.11 \
                Supercontig_12.12 Supercontig_12.13 Supercontig_12.14 Supercontig_12.15 Supercontig_12.16 Supercontig_12.17 \
                Supercontig_12.18 Supercontig_12.19 Supercontig_12.20 \
            --group-file Cross_E/Cross_E_s40.bam_groups.txt \
            --min-group-size 1 --max-group-size 1 | grep -v "#" \
            >> Cross_E/Cross_E_s40.simulated.dat

<br />



#### Step 3: extract synthetic read pairs

        ## Extract synthetic nucleotide changes
        awk 'BEGIN{OFS="\t"} /^\#Tag/ || $1 == "MUT" {if(/\#Tag/){$2 = "#Chrom";} print;}' \
            Cross_E/Cross_E_s40.simulated.dat | \
            cut -f 2- | sort -k1,1 -k2,2n \
            > Cross_E/Cross_E_s40.simulated.vars.csv
        
        ## prepare a vcf-format file
        awk 'BEGIN{OFS="\t";} {if(/#Chrom/) {
                print "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tMUTATION";next;}
                print $1,$2,$4,$5,$6,"999\t.\tTYPE=snp\tSDP",$7;}' \
            Cross_E/Cross_E_s40.simulated.vars.csv | \
            cat NC.fileheader.vcf - \
            > Cross_E/Cross_E_s40.simulated.vars.vcf
        
        ## check how many synthetic mutations in duplicates/duplicate-proximal/non-duplicates
        awk 'BEGIN{OFS="\t"} {if(!/\#/){$10 = ".";} print;}' \
            Cross_E/Cross_E_s40.simulated.vars.vcf | \
            bcftools annotate -a duplicate_regions.tab.gz \
            -h repeat_status.hdr -c CHROM,FROM,TO,RepeatState - \
            > Cross_E/Cross_E_s40.simulated.vars.regions.vcf
            
        
        ## generate a bed-format file
        awk 'BEGIN{OFS="\t";} !/\#/ {print $1,($2-1),$2;}' \
            Cross_E/Cross_E_s40.simulated.vars.vcf \
            > Cross_E/Cross_E_s40.simulated.vars.bed
        
        ## Extract reads carry the synthetic mutations
        awk '$1 ~ /SAM:/' Cross_E/Cross_E_s40.simulated.dat | \
            cut -f 2- \
            > Cross_E/Cross_E_s40.simulated.reads.sam


        ## Generate synthetic reads
        find  -name "*.bam" -print | xargs -n 1 -P 4 -I BAM_FILE sh -c '
                sample=`basename BAM_FILE | cut -d"." -f1`
                
                echo "${sample}"
                
                extract_bam_pairs.pl --samtools "-F 3844" --extend 2000 --bam BAM_FILE \
                    --input Cross_E/Cross_E_s40.simulated.vars.csv \
                    --patches Cross_E/Cross_E_s40.simulated.reads.sam \
                    > Cross_E/reads/${sample}_ex2k.sam \
                    2> Cross_E/reads/${sample}_ex2k.log
                
                samtools view -H BAM_FILE | \
                    cat - Cross_E/reads/${sample}_ex2k.sam \
                    > Cross_E/reads/${sample}_ex2k.add.sam
                
                ## generate a bam file of extracted reads for debugging purpose
                java -jar picard-tools-1.114/SortSam.jar \
                    INPUT=Cross_E/reads/${sample}_ex2k.add.sam \
                    OUTPUT=Cross_E/reads/${sample}_ex2k.sort.bam \
                    SORT_ORDER=coordinate VALIDATION_STRINGENCY=LENIENT \
                    > Cross_E/reads/${sample}_ex2k.sort.log 2>&1
                
                samtools index Cross_E/reads/${sample}_ex2k.sort.bam
                
                ## convert to fastq file
                sam2fastq.pl -i Cross_E/reads/${sample}_ex2k.sam \
                    -o Cross_E/reads/${sample}_ex2k \
                    >> Cross_E/reads/${sample}_ex2k.log 2>&1 && \
                    rm -v Cross_E/reads/${sample}_ex2k.sam \
                          Cross_E/reads/${sample}_ex2k.add.sam
            '
        
        find Cross_E/reads -name "*.fq" | \
            xargs -n 1 -P 6 -I {} gzip -v {}


        ## check extracted read pairs and make sure the number of extracted read pairs are
        ## roughly in proportion with the original sequencing depth, otherwise try to re-run
        ## the extraction steps
        grep "pairs were found" Cross_E/reads/*_ex2k.log | \
            awk '{print $1"\t"$2"\t"$10;}' \
            > Cross_E/Cross_E_s40.simulated.read_pairs.log

<br />


#### Step 4: remap the simulated read pairs to reference genome

        find Cross_E/reads/ -name "*_1.fq.gz" | \
            sed 's/_1.fq.gz$//' | xargs -n 1 -P 4 -I PREFIX \
            sh -c '
                sample=`basename PREFIX | sed "s/_ex2k//"`
                lane_id=`basename PREFIX | sed "s/_ex2k//"`
                
                echo "[`date`]: Start mapping ${sample}:${lane_id} ... "
                
                read1=PREFIX"_1.fq.gz"
                read2=PREFIX"_2.fq.gz"
                
                ## Align reads with BWA-MEM algorithm
                bwa mem -t 1 -M -R "@RG\tID:${lane_id}\tLB:${sample}\tPL:Illumina\tPU:${sample}\tSM:${sample}" \
                    reference_genome/neurospora_crassa.fasta ${read1} ${read2} \
                    > Cross_E/mapped/${sample}.NC12.bwa.sam \
                    2> Cross_E/mapped/${sample}.NC12.bwa.process.log
                
                ## sort bam file
                java -jar picard-tools-1.114/SortSam.jar \
                    INPUT=Cross_E/mapped/${sample}.NC12.bwa.sam \
                    OUTPUT=Cross_E/mapped/${sample}.NC12.bwa.sort.bam \
                    SORT_ORDER=coordinate VALIDATION_STRINGENCY=LENIENT \
                    >> Cross_E/mapped/${sample}.NC12.bwa.process.log 2>&1 && \
                    rm -v Cross_E/mapped/${sample}.NC12.bwa.sam
                
                echo "[`date`]: Start marking duplicates ${sample} ... "
                
                ## mark duplicates
                java -jar picard-tools-1.114/MarkDuplicates.jar \
                    MAX_FILE_HANDLES_FOR_READ_ENDS_MAP=8000 VALIDATION_STRINGENCY=LENIENT \
                    INPUT=Cross_E/mapped/${sample}.NC12.bwa.sort.bam \
                    OUTPUT=Cross_E/mapped/${sample}.NC12.bwa.sort.dedup.bam \
                    METRICS_FILE=Cross_E/mapped/${sample}.NC12.bwa.sort.dedup.metrics \
                    >> Cross_E/mapped/${sample}.NC12.bwa.process.log 2>&1 && \
                    rm -v Cross_E/mapped/${sample}.NC12.bwa.sort.bam
                
                ## index bam file
                samtools index Cross_E/mapped/${sample}.NC12.bwa.sort.dedup.bam
                
                
                echo "[`date`]: Start realigning ${sample} ... "
                
                ## realignment
                java -jar GenomeAnalysisTK-3.7/GenomeAnalysisTK.jar \
                    -R reference_genome/neurospora_crassa.fasta \
                    -T RealignerTargetCreator -nt 1 \
                    -o Cross_E/mapped/${sample}.NC12.bwa.sort.dedup.realn.intervals \
                    -I Cross_E/mapped/${sample}.NC12.bwa.sort.dedup.bam \
                    >> Cross_E/mapped/${sample}.NC12.bwa.process.log 2>&1
                
                java -jar GenomeAnalysisTK-3.7/GenomeAnalysisTK.jar \
                    -R reference_genome/neurospora_crassa.fasta \
                    -T IndelRealigner \
                    -targetIntervals Cross_E/mapped/${sample}.NC12.bwa.sort.dedup.realn.intervals \
                    -o Cross_E/mapped/${sample}.NC12.bwa.sort.dedup.realn.bam \
                    -I Cross_E/mapped/${sample}.NC12.bwa.sort.dedup.bam \
                    >> Cross_E/mapped/${sample}.NC12.bwa.process.log 2>&1 && \
                    rm -v Cross_E/mapped/${sample}.NC12.bwa.sort.dedup.bam \
                          Cross_E/mapped/${sample}.NC12.bwa.sort.dedup.bam.bai
                
                echo "[`date`]: Finished processing ${sample}"
            '

<br />

#### Step 5: Variants discovery

            
        ## HaplotypeCaller joint calling across the genome
        BAM_FILEs=`find Cross_E/mapped/ -name "*.bam" -print | \
            xargs -I BAM_FILE echo -n "-I BAM_FILE "`
        java -jar GenomeAnalysisTK-3.7/GenomeAnalysisTK.jar \
            -R reference_genome/neurospora_crassa.fasta \
            -T HaplotypeCaller -stand_call_conf 30.0 \
            -o Cross_E/variants/Cross_E_s40_sim.hc.genome.vcf \
            ${BAM_FILEs} 2>&1 | \
            tee Cross_E/variants/Cross_E_s40_sim.hc.genome.log
        
        
        ## HaplotypeCaller joint calling in target regions
        ## This is to avoid unusual behavior of HaplotypeCaller when only part of reads given, HaplotypeCaller occasionally lost some variants due to unknown reasons
        BAM_FILEs=`find Cross_E/mapped/ -name "*.bam" -print | \
            xargs -I BAM_FILE echo -n "-I BAM_FILE "`
        java -jar GenomeAnalysisTK-3.7/GenomeAnalysisTK.jar \
            -R reference_genome/neurospora_crassa.fasta \
            -T HaplotypeCaller -stand_call_conf 30.0 \
            -L Cross_E/Cross_E_s40.simulated.vars.bed \
            -o Cross_E/variants/Cross_E_s40_sim.hc.sim_only.vcf \
            ${BAM_FILEs} 2>&1 | \
            tee Cross_E/variants/Cross_E_s40_sim.hc.sim_only.log
        
        
        ## combine both results
        vcf_process.pl --vcf Cross_E/variants/Cross_E_s40_sim.hc.genome.vcf \
            --secondary-vcf Cross_E/variants/Cross_E_s40_sim.hc.sim_only.vcf \
            --combine-rows 0 1 --compare-rows 2 --primary-tag Genome --secondary-tag Target --intersect-tag "Genome+Target" \
            > Cross_E/variants/Cross_E_s40_sim.hc.combined.vcf
        
        
        ## check how many sites are callable
        map_records.pl --subject Cross_E/variants/Cross_E_s40_sim.hc.combined.vcf \
            --query Cross_E/Cross_E_s40.simulated.vars.csv \
            --rows1 0 1 --rows2 0 1 | cut -f 1,2,11,14 | \
            perl -ne 'if (/\#/){print; next;} my @line = (split /\s+/);
                unless($line[2]) {$line[2] = "-"; $line[3] = "-";}; $line[3] =~ s/.*Combine=//;
                my $out_line = join "\t", @line; print "$out_line\n";' \
            > Cross_E/variants/Cross_E_s40_sim.hc.combined.called.csv
        
        
        ## check pre-existing variants, a simulated mutation same as pre-existing variant is invalid as it could never be recovered
        bcftools view -s 2489,4200 \
            --min-ac 1 Cross_E/variants/Cross_E_s40_sim.hc.combined.vcf \
            > Cross_E/variants/Cross_E_s40_sim.hc.controls.vcf
        
        map_records.pl --subject Cross_E/variants/Cross_E_s40_sim.hc.controls.vcf \
            --query Cross_E/Cross_E_s40.simulated.vars.csv \
            --rows1 0 1 --rows2 0 1 | \
            perl -ne 'if(/\#Chrom/){print "#CHROM\tPOS\tIs_Exist\tVAR(Control)\n";} next if (/\#/);
                my @line = (split /\s+/); my $tag = "NO";
                if($line[8] ne "N/A") {if($line[10] eq $line[5]){$tag = "YES";}else {$tag = "DIFF";}} else {$line[10] = "-";}
                print "$line[0]\t$line[1]\t$tag\t$line[10]\n";' \
            > Cross_E/variants/Cross_E_s40_sim.hc.controls.called.csv
            
            
<br />

#### Step 6: Mutation detection
            
            
        detect_tetrad_mutations.pl -v Cross_E/variants/Cross_E_s40_sim.hc.combined.vcf \
            -g Cross_E/Cross_E_s42.sample_groups.txt \
            --min-supp-depth 5 \
            > Cross_E/hc_gvcf/Cross_E_s40_sim.hc.combined.mut.vcf
        
        
        ## compare recovered results with original synthetic sites
        map_records.pl --subject Cross_E/hc_gvcf/Cross_E_s40_sim.hc.combined.mut.vcf \
            --query Cross_E/Cross_E_s40.simulated.vars.csv \
            --rows1 0 1 --rows2 0 1 | cut -f 3,10 --complement | \
            perl -ne 'if (/\#/){print; next;} my @line = (split /\s+/); $line[11] =~ s/.*Combine=//;
                if ($line[7] eq "N/A") {$line[7] = "-"; $line[8] = "-"; $line[9] = "-";  $line[11] = "-"; $line[12] = "-"; $line[13] = "-";}
                if ($line[2] =~ /\;/){$line[10] = "Sexual-2:2";} else {$line[10] = "Sexual-3:1";}
                my $out_line = join "\t", @line; print "$out_line\n";' \
            > Cross_E/hc_gvcf/Cross_E_s40_sim.hc.combined.mut.called.csv
            
            


<br />




