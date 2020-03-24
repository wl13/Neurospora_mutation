# input dir   f19-mrout
opendir(DIR,"$ARGV[0]");
foreach my $file (readdir DIR) {
	if ($file=~/gz/) {
		open FILE, "gzip -dc $ARGV[0]\/$file | " or die "Can not open file $file\n";
		while (<FILE>) { s/\r//g;  chomp;
			($chr1,$start1,$end1,$chr2,$start2,$end2) = split /\D+/;
			if ($chr1 == $chr2 and $start1 != $start2)  { 
				$hash{"$chr1\-$chr2\:$start1\-$start2"} = "$chr1\-$chr2\:$end1\-$end2"; 
			}
#			if (($chr1 != $chr2) or ($chr1 == $chr2 and $start1 != $start2)) { $hash{"$chr1\-$chr2\:$start1\-$start2"} = "$chr1\-$chr2\:$end1\-$end2"; }
				#10:100185..100238	10:100185..100240
				#10:100185..100238	2:3312562..3312617
		}
		close FILE;
		$line ++; $time = localtime(); warn "$time processed $line files $file in stage 1\n";
	}
}
# foreach $t (keys %hash) {print "$t\t$hash{$t}\n";}
$line=0;
open G11, ">g11.homochr-out";
open G16, ">g16.homochr-out";
open G21, ">g21.homochr-out";

opendir(DIR,"$ARGV[0]");
foreach my $file (readdir DIR) {
	if ($file=~/gz/) {
		open FILE, "gzip -dc $ARGV[0]\/$file | " or die "Can not open file $file\n";
		while (<FILE>) { s/\r//g;  chomp;
			($chr1,$start1,$end1,$chr2,$start2,$end2) = split /\D+/;
	#		for ($i=0;$i<=$#c-1;$i++) {
	#			for ($j=$i+1;$j<=$#c;$j++) {
				if ($chr1 == $chr2 and $start1 != $start2)  {
					$t = "$chr1\-$chr2\:$start1\-$start2"; # push @new, $new;
					if (exists $hash{$t}) { 
						$g11 = $hash{$t}; 
						if (exists $hash{$g11}) {
							$g16 = $hash{$g11};  print G11 "$t\t$g11\t$g16\n"; 
							if (exists $hash{$g16}) {
								$g21 = $hash{$g16}; print G16 "$t\t$g11\t$g16\t$g21\n";
								if (exists $hash{$g21}) {
									print G21 "$t\t$g11\t$g16\t$g21\t$hash{$g21}\n";  
								}  
							}
						}
					}
				}
	#		}
		}
		$line ++; $time = localtime(); warn "$time processed $line files $file in stage 2\n";
		close FILE;
	}
}
close G11;  close G16; close G21;
__DATA__
foreach $t (sort keys %hash) {
	if (exists $hash{$t}) { $call{$t} = $hash{$t};}
}
foreach $t (sort keys %call){
	print "$t\t$call{$t}\n";
}
