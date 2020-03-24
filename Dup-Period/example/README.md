
#### Step1: chopping genome into short fragments/units with all combinations of 10~12bp
* Read into the mock_genome and directly chop into short fragments, write to folder "f08-u3g6"
	
		mkdir f08-u3g6
		
		perl ./1-1.chop_genome_to_group-ff.pl


#### Step2: searching for units with at least 3 identical base in front
* The results will write to folder "f12"

		mkdir f12
		
		find ./f08-u3g6 -maxdepth 1 -name "g-*.fas.gz" | \
		    sed 's/.fas.gz$//' | sed 's/g-//' |  xargs -n 1 -P 2 -I PREFIX \
		    sh -c '
		        sample=`basename PREFIX`
			
		    	echo "[`date`]: Start searching ${sample} ... "
				perl ./2-1.search_matched_units-ff.pl ${sample}
		    '


#### Step3 : linking matched units in same chromosomes
* This step generates three files g11.homochr-out, g16.homochr-out, and g21.homochr-out, these are interspersed repeats identified in given genome
	
		perl ./3-1.join_matched_units_in_same_chr.pl f12


#### Step4: gathering joined units if they satisfy the specied criterion like >= 16 units, and write to bed output
* The final results are given in bed format after merging all interspersed repeat regions, the right boundary is extended to 50bp downstream so might exceed the chromosome length if reach the end
	
		cat ./g11.homochr-out \
		    ./g16.homochr-out \
		    ./g21.homochr-out | \
		    perl ./4.gather_joined_units.pl \
		    > final.bed


