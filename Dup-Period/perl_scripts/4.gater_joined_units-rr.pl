open LEN, "</mnt/san1/usr/thw/info/neurospora/chr-length.txt" or die "Can not open genome file\.\n";
while (<LEN>) { s/\r//g; chomp; 
	@c = split /\t/;
	$len{$c[0]} = $c[1];
}

while (<>) { s/\r//g; chomp;
	@c = split /\t/;	($chr0,$chr1,$start0,$start1) = split /\D/, $c[0]; ($chr0,$chr1,$end0,$end1) = split /\D/, $c[-1];
	$chr0="Supercontig\_12\.$chr0"; $chr1="Supercontig\_12\.$chr1"; 

	# active the next two lines if active reverse fasta seq
	$start1 = $len{$chr1} - $start1; $end1 = $len{$chr1} - $end1;
	$tmp = $start1; $start1 = $end1; $end1 = $tmp;

	$cha = $start1 - $start0; $cha = 0 - $cha if $cha < 0;
	if ( $chr0 ne $chr1 or   (($chr0 eq $chr1) and $cha >=150 ) ) {
		for ($i=$start0; $i<=$end0+50; $i++) { $hash{$chr0}{$i} ++; }
		for ($i=$start1; $i<=$end1+50; $i++) { $hash{$chr1}{$i} ++; }
	}
}

# $ss = keys %count;

foreach $chr ( sort keys %hash ) {
#	print "xxx $t\n";
	my $hash2 = $hash{$chr}; $keycount = keys %$hash2; warn "keycount $keycount\t";
	$start = 0; $stop = 0; $startex = 0; $stopex = 0; $i = 0;
	foreach $x ( sort { $a <=> $b } keys %$hash2 ) {
#		print "$x\t".$hash2->{${x}+1}."\n";
		if ( exists $hash2->{${x}-1} ) {
			$stop = $x ;			
		} else {
			unless ($start == 0) {  $length = $stop - $start + 1; print "$chr\t$start\t$stop\tdupT\t$length\n"; $accu += $length; }
#			print "$chr\t".$hash2->{$x}."\tstart $start\tstop $stop\n";
			$start = $x;
		}
		$i ++; if ( $i == $keycount ) { $length = $stop - $start + 1; warn "yes\n"; print "$chr\t$start\t$stop\tdupT\t$length\n"; $accu += $length; }
	}
}
warn "$single\t$ss\t$accu\n";
__DATA__
foreach $t ( sort keys %hash ) {
	print "xxx $t\n";
	my $hash2 = $hash{$t};
	foreach $x ( sort { $a <=> $b } keys %$hash2 ) {
		print "$x\t".$hash2->{$x}."\n";
	}
}
