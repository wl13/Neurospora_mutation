open GENOME, "< mock_genome.rename.fasta" or die "can not open genome file";
while (<GENOME>) { s/\r//g; chomp;
	if (/>(.+)/){
		$id = (split /\s+/, $1)[0];
	} else {
		$fas{$id} .= $_;
	}
}
###foreach $t ( keys %fas ) {
###	$fas{$id} = reverse $fas{$id};
###	$fas{$id} =~tr/ATCG/TAGC/;
###}

open LEN, "<chr-length.txt" or die "Can not open genome file\.\n";
while (<LEN>) { s/\r//g; chomp; 
	@c = split /\t/;
	$len{$c[0]} = $c[1];
}
open ARRAY, "<unit_length_combinations.txt" or die "Can not open array file\.\n";
while (<ARRAY>) { chomp;
	@c = split /\t/; $group = join "\_", @c;    #$group1 = $c[0]; $i++; $iout = sprintf "%04d", $i;
	foreach $chr ( sort keys %fas ) {
		$seq = $fas{$chr};	$maxlen = $len{$chr};
#		warn "Initail task for chromosome $chr containing $maxlen bp\.\n";
		for ( $segstart=0; $segstart<= $maxlen-150; $segstart++ ) {
			$start=$segstart; 
			$new = substr($seq,$start,3);
			foreach $c (@c) {	$start+=$c; $new .= substr($seq,$start,3);  }
			$newn = substr($new,0,5);
			$hash{$newn} .= ">$chr\-$segstart\-$group\n$new\n" unless $new =~/N/;
		}
	}
	$ij++;  $time = localtime(); warn "$time processed array $ij $group\n";
}
foreach $t ( keys %hash ) {
	open GROUP, " | gzip > f08\-u3g6\/g\-$t\.fas\.gz" or die $!;
	print GROUP $hash{$t};
	close GROUP;
	delete $hash{$t};
}
