
### Scripts and pipelines for mutation analyses in Neurspora crassa

<br />

* An example from preprocessing to mutation calling could be found in Simulation Step 4~6

## Scripts

#### detect_tetrad_mutations.pl   
> Detect mutations in parent-tetrad/ascus samples

* **Usage:**

        detect_tetrad_mutations.pl -v variants.vcf -g sample_groups.txt --min-supp-depth 5 mut.vcf


* **Options:**   

        -v, --vcf     <filename>
            input vcf file, required
            
        Note: This script is designed to process vcf files with AD (Allele Depth) field for each sample
         
        -o, --output  <filename>
            output filename, default to STDOUT
        
        -q, --quality     <float>
            loci with quality smaller than this value will filtered
        
        
        -g, --group-file  <file>
            file contain group infos of each sample, each sample per line, e.g.
            sample1 group1
            sample2 group1
            sample3 group2
            ...
            set this option to screen group-specific mutation, only samples belong
            to different groups would be used as compare samples
            
        
        --min-supp-depth  <int>
            minimum number of supporting reads, for multiple samples carry the same
            mutation, only one need to pass this criterion [default: 1]
    

<br />


