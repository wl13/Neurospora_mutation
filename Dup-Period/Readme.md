
### Code to search for "Dup-Period" regions, i.e., regions with weak interspersed homology which follow the matching period defined by Gladyshev et. al. [ref 1,2].


* A matching period of 10~12 with at least 3 matching bases in the beginning of each period (e.g., 3H7N, 3H8N, 4H7N, etc.) were used. A requirement of 16 such matching units are required here (i.e. over 160 bp in length), though decreasing (till ~10) or increasing (till ~21) this threshold differs few in final results. 

* Original code contributed by Huawei Tan


* References

    1) Gladyshev E, Kleckner N. Direct recognition of homology between double helices of DNA in Neurospora crassa. Nature Communications. 2014;5:3509.
    2) Gladyshev E, Kleckner N. Recombination-Independent Recognition of DNA Homology for Repeat-Induced Point Mutation (RIP) Is Modulated by the Underlying Nucleotide Sequence. PLOS Genet. 2016;12:e1006015. 


<br />


* The mentioned perl scripts are avaliable in the folder "perl_scripts". These example scripts are given for a single strand and orientation, one could easily modify them for other strand or orientation. Note the current scripts need much improvement for generally use as well as better efficiency. A full example on a mock genome could be found in the "example" folder

#### Step1 (1.chop_genome_to_group-rr.pl): chopping genome into short fragments/units with all combinations of 10~12bp

* All 243 (=3^5) possible combinations of lengths are given in file "unit_length_combinations.txt"


#### Step2 (2.search_matched_units.pl): searching for units with at least 3 identical base in front 


#### Step3 (3-1.join_matched_units_in_same_chr.pl): linking matched units in same or different (3-2.join_matched_units_in_diff_chr.pl) chromosomes


#### Step4 (4.gather_joined_units.pl): gathering joined units if they satisfy the specied criterion like >= 16 units, and write to bed output



<br />


