$dir1 = 'f16-minusrev-u3g6';
$dir2 = 'f08-u3g6';
$dirout = 'f19-mrout';

$sample=$ARGV[0];
#warn "$sample";
	
open FILE, "gzip -dc $dir1\/mr\-${sample}\.fas\.gz | " or die "Can not open seq file\n";
while (<FILE>) { s/\r//g; chomp;
	if (/>(.+)/) {
		($chr1,$start1,$period1) = split /\-/, $1; $chr1=~s/Supercontig\_12\.//;
		@period1 = split /\_/, $period1; $end1 = $start1; foreach $p1 (@period1) { $end1 += $p1; }
		$id1 = "$chr1\:$start1\.\.$end1";
	} else {
		$hash{$_} .= "$id1\t";
	}
}
close FILE;

open FILE, "<$dir2\/g\-$sample\.fas" or die "Can not open seq file\n";
open FILEOUT, " | gzip > $dirout\/mrp\-$sample\.txt\.gz" or die "Can not open seq file\n";
while (<FILE>) { s/\r//g; chomp;
	if (/>(.+)/) {
		($chr2,$start2,$period2) = split /\-/, $1; $chr2=~s/Supercontig\_12\.//;
		@period2 = split /\_/, $period2; $end2 = $start2; foreach $p2 (@period2) { $end2+= $p2; }
		$id2 = "$chr2\:$start2\.\.$end2";
	} else {
		$seq=$_;
		if ( exists $hash{$seq} ) { $name2 = $hash{$seq}; chop $name2; @name2 = split /\t/, $name2;
			foreach $nn (sort @name2) { $new{"$id2\t$nn"} ++ if $nn ne $id2;}
			#print FILEOUT "$id2\t$nn\n"; 
		}
	}
}
close FILE;

foreach $t (sort keys %new) { print FILEOUT "$t\n"; }
close FILEOUT;



__DATA__
opendir(DIR, "$ARGV[0]") or die "Can not open dir \n";
foreach my $file (readdir DIR) {
if ($file=~/fas/) {
	undef %hash; #warn $file;
	$fileout = (split /\W/, $file)[1];
	$line ++; if ($line%10 ==0) { $time = localtime(); warn "$time start to process $line files $file \n"; }
	open FILE, "<$ARGV[0]\/$file" or die "Can not open seq file\n";
	while (<FILE>) { s/\r//g; chomp;
		if (/>(.+)/) {
			$name = $1;
		} else {
			$seq=$_; $seq=tr/ATCG/TAGC/; $hash{$name} = $seq;
		}
	}
	close FILE;
	open FILE, "<$ARGV[0]\/$file" or die "Can not open seq file\n";
	open FILEOUT, ">$ARGV[1]\/minusfor\-$fileout\.txt" or die "Can not open seq file\n";
	while (<FILE>) { s/\r//g; chomp;
		if (/>(.+)/) {
			($chr,$start,$period) = split /\-/, $1; $chr=~s/Supercontig\_12\.//;
			@period = split /\_/, $period; $end = $start; foreach $p (@period) { $end+= $p; }
			$id = "$chr\:$start\.\.$end";
		} else {
			$seq=$_; 
			foreach $t ( sort keys %hash) {
				if ( $hash{$t} eq $seq) {
					($chr2,$start2,$period2) = split /\-/, $t; $chr2=~s/Supercontig\_12\.//;
					@period2 = split /\_/, $period2; $end2 = $start2; foreach $p2 (@period2) { $end2+= $p2; }
					$id2 = "$chr2\:$start2\.\.$end2";
					print FILEOUT "$id\t$id2\n";
				}
			}
		}
	}
	close FILE;  close FILEOUT;
}
}

#		$cha=$start1-$start0; $cha = 0 - $cha if $cha <0;
#		if ( $chr0 ne $chr1 or   (($chr0 eq $chr1) and $cha >=150 ) ) { chop $hash2{$t}; print "$hash2{$t}\n"; }


