#!/usr/bin/perl -w
#
#   detect_mutations.pl -- Screen out candidate mutation sites.
#                          
#
#   Author: Nowind
#   Created: 2012-05-31
#   Updated: 2019-07-12
#   Version: 1.0.0
#
#   Change logs:
#   Version 1.0.0 19/07/12: The initial version.


use strict;

use Data::Dumper;
use Getopt::Long;
use File::Find::Rule;
use File::Basename;
use List::Util qw( sum max );

use MyPerl::FileIO qw(:all);
use MyPerl::Vcf qw(:all);

######################## Main ########################
my $CMDLINE = "perl $0 @ARGV";
my $VERSION = '1.0.0';
my $HEADER  = "##$CMDLINE\n##Version: $VERSION\n";
my $SOURCE  = (scalar localtime()) . " Version: $VERSION";

my %options = ();
   $options{max_cmp_depth}  = 1;
   $options{help} = ($CMDLINE =~ /\-help/) ? 0 : 1;
GetOptions(
            "vcf=s"              => \$options{vcf},
            "output=s"           => \$options{output},
            
            "qual=f"             => \$options{min_qual},
            
            "min-supp-depth=i"   => \$options{min_supp_depth},
            "min-supp-plus=i"    => \$options{min_supp_plus},
            "min-supp-minus=i"   => \$options{min_supp_minus},
            
            "min-lib-cnt=i"      => \$options{min_lib_cnt},
            "min-lib-depth=i"    => \$options{min_lib_depth},
            
            "max-cmp-miss=i"     => \$options{max_cmp_missing},
            "max-cmp-depth=i"    => \$options{max_cmp_depth},
            "max-cmp-perc=f"     => \$options{max_cmp_perc},
            "max-cmp-total=i"    => \$options{max_cmp_total},
            "max-cmp-freq=s"     => \$options{max_cmp_freq},
            
            "max-shared-freq=i"  => \$options{max_shared_freq},
            "group-file=s"       => \$options{group_file},
            
            
            "no-ref-mut"         => \$options{no_ref_mut},
            
            "mask-only=s{,}"     => \@{$options{mask_only}},
            "controls=s{,}"      => \@{$options{control_samples}},
            
            "min-indel-len=i"    => \$options{min_indel_len},
            "max-indel-len=i"    => \$options{max_indel_len},
            
            "append-info"        => \$options{append_info},
           );


unless( $options{vcf} && $options{help} ) {
    print <<EOF;

$0  -- Screen out candidate mutation sites.

Version: $VERSION

Usage:   perl $0 [options]

Options:
    -v, --vcf     <filename>
        input vcf file, required
        
    *Note: This script is designed to process vcf files with AD (Allele
     Depth) field for each sample
     
    -o, --output  <filename>
        output filename, default to STDOUT
    
    -a, --append-info
        Append new INFO field without replacing the original one 
    
    -f, --filter  <strings>
        skip filter loci, can have multiple values, separate by space, e.g.
        "LowQual SNPFilter ..."
    -M, --match   <strings>
        only retain loci matches, can have multiple values, separate by space,
        e.g. "PASS ..."


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
        
    --max-shared-freq <int>
        locus with allele frequency below this value will be considering as a
        shared mutation locus
    
    --min-supp-depth  <int>
        minimum number of supporting reads, for multiple samples carry the same
        mutation, only one need to pass this criterion [default: 1]

    
    --controls  <strings>
        specify samples served as controls where no missing calls is allowed,
        and shared mutations contain those samples will be filtered

EOF

    exit(1);
}

$|++;



if ($options{output}) {
    open (STDOUT, "> $options{output}") || die $!;
}


print STDERR "# $0 v$VERSION\n# " . (scalar localtime()) . "\n";



print STDERR ">> Start detecting candidate mutations in $options{vcf} ... ";
detect_mutations(\%options);
print STDERR "done!\n";


print STDERR "# " . (scalar localtime()) . "\n";

######################### Sub #########################

=head2 get_group_info

    About   : Get group infos of each sample
    Usage   : get_group_info($group_file);
    Args    : File contain group infos
    Returns : Null

=cut
sub get_group_info
{
    my ($in, $rh_group_infos) = @_;
    
    my $fh = getInputFilehandle($in);
    while (<$fh>)
    {
        next if (/\#/ || /^\s+$/);
        
        my ($sample_id, $group_id) = (split /\s+/);
        
        $rh_group_infos->{$sample_id} = $group_id;
    }
}

=head2 check_genotype

    About   : Check homozygous/heterozygous/missing by GT
    Usage   : check_genotype($GT);
    Args    : GT info
    Returns : Null

=cut
sub check_genotype
{
    my ($GT) = @_;
    
    if ($GT !~ /\||\//) {
        return "MISS";
    }
    
    my ($allele1, $allele2) = (split /\||\//, $GT);
    
    if ($allele1 eq ".") {
        return "MISS";
    }
    elsif ($allele1 == $allele2) {
        return "HOM";
    }
    else {
        return "HET";
    }
}


=head2 detect_mutations

    About   : Detect candidate mutations
    Usage   : detect_mutations($vcf_file);
    Args    : Vcf file contains all samples
    Returns : Null

=cut
sub detect_mutations
{
    my ($opts) = @_;
    

    ##
    ## parse group infos
    ##
    my %group_infos = ();
    if ($opts->{group_file}) {
        ###print STDERR "#DEBUG: Reading $opts->{group_file}#\n";
        get_group_info($options{group_file}, \%group_infos);
    }

    
    ##
    ## set default values
    ##
    $opts->{min_supp_depth}  ||= 0;
    
   
    my @sample_ids  = ();
    my %sample_rows = ();
    my $out_header  = '';
    my $fh = getInputFilehandle($opts->{vcf});
    while (<$fh>)
    {
        if (/#CHROM/) {
            my @vcf_header = (split /\s+/);
            
            next if (@sample_ids > 0);
            
            @sample_ids = @vcf_header[9..$#vcf_header];
            
            for (my $i=0; $i<@sample_ids; $i++)
            {
                $sample_rows{$sample_ids[$i]} = $i;
            }
            
            $out_header .= <<EOF;
##INFO=<ID=MA,Number=1,Type=String,Description="Mutation allele">
##INFO=<ID=MAR,Number=1,Type=Float,Description="Ratio of reads contain mutation allele among all covered reads in mutation sample">
##INFO=<ID=FPD,Number=1,Type=Integer,Description="Depth of mutation-like allele in compare samples">
##INFO=<ID=FPFQ,Number=1,Type=Integer,Description="Mutation-like allele frequency of different depth in compare samples">
##INFO=<ID=FPS,Number=1,Type=String,Description="Compare samples with mutation-like alleles">
##INFO=<ID=GRPID,Number=1,Type=Integer,Description="Group ID">
##INFO=<ID=GRPD,Number=1,Type=Integer,Description="Depth of mutation-like allele in other group members">
##INFO=<ID=GRPFQ,Number=1,Type=Integer,Description="Mutation-like allele frequency of different depth in other group samples">
##INFO=<ID=GRPS,Number=1,Type=String,Description="Other group samples with mutation-like alleles">
##INFO=<ID=NMISS,Number=1,Type=Integer,Description="Number of uncallable compare samples">
##INFO=<ID=SMISS,Number=1,Type=String,Description="Uncallable compare samples">
##INFO=<ID=Shared,Number=.,Type=String,Description="Number of samples sharing this mutation allele(Details of shared samples, listed same as FORMAT field)">
##source=$SOURCE $CMDLINE
EOF
            print "$out_header";
            print "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tMUTATION\n";
            
            next;
        }
        elsif (/\#\#/ || /^\s+$/) {
            $out_header .= $_; next;
        }
        
        my ($CHROM, $POS, $ID, $REF, $ALT, $QUAL, $FILTER,
            $INFO, $FORMAT, @samples_all) = (split /\s+/);
        
        ###next unless($POS == 84943);
        
        next if (defined $opts->{min_qual} && ($QUAL eq '.' || $QUAL < $opts->{min_qual}));   ## filter by Quality
        
        my @vars = ($REF, (split /\,/, $ALT));
        my @tags = (split /\:/, $FORMAT);
        my %tags = ();
        for (my $i=0; $i<@tags; $i++) { $tags{$tags[$i]} = $i; }
        
        unless(defined($tags{AD})) {
            print STDERR "Warning: no AD tag found in line $_, skip\n"; next;
        }
        
        ##
        ## parse read depth
        ##
        my %read_counts     = ();
        my %allele_freq     = ();
        my %sample_infos    = ();
        for (my $i=0; $i<@samples_all; $i++)
        {
            my ($GT, $AD) = (split /\:/, $samples_all[$i])[$tags{GT}, $tags{AD}];
            
            if (!$AD || ($AD eq ".") || ($AD eq "./.")) {
                $AD = 0;
            }
            
            my $group_id = $group_infos{$sample_ids[$i]};
            push @{$sample_infos{$group_id}->{GT}}, $GT;
            push @{$sample_infos{$group_id}->{DP}}, sum((split /\,/, $AD));
            
            push @{$sample_infos{$group_id}->{index}}, $i;
        }
        
        ## make sure no control sample is heterozygous
        my $is_control_het     = grep { check_genotype($_) eq "HET"; } @{$sample_infos{control}->{GT}};
        
        next if ($is_control_het > 0);
        
        
        ## make sure at least one control sample is properly covered
        my $is_control_lowdepth = grep { $_ >= $opts->{min_supp_depth}; } @{$sample_infos{control}->{DP}};
        
        next if($is_control_lowdepth < 1);
        
        ###print STDERR Dumper(%sample_infos); exit;
        
        my @mutated_samples = ();
        my @mutated_infos   = ();
        
        for my $group_id (sort keys %sample_infos)
        {
            next if ($group_id eq "control");
            
            my @spores_diff = ();
            my $spores_pass = 0;
            
            my @tetrad_spores = @{$sample_infos{$group_id}->{GT}};
            
            for (my $j=0; $j<@tetrad_spores; $j++)
            {
                my $is_same_as_control = grep { $_ eq $tetrad_spores[$j] } @{$sample_infos{control}->{GT}};
                
                if ($is_same_as_control == 0) {
                    push @spores_diff, $sample_infos{$group_id}->{index}->[$j];
                    
                    if ($sample_infos{$group_id}->{DP}->[$j] >= $opts->{min_supp_depth}) {
                        $spores_pass ++;
                    }
                }
            }
            
            if ($spores_pass > 0) {
                push @mutated_samples, @sample_ids[@spores_diff];
                push @mutated_infos,   @samples_all[@spores_diff];
            }
        }
            
        next unless(@mutated_samples > 0);
        
        my $mut_samples  = join ";", @mutated_samples;
        my $mut_infos    = join ":", @mutated_infos;
        
        print "$CHROM\t$POS\t$mut_samples\t$REF\t$ALT\t$QUAL\t$FILTER\t$INFO\t$FORMAT\t$mut_infos\n";
    }
}

