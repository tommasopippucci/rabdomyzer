#!/usr/bin/env perl
                                                                                                   # < When family, no gene. When gene, no family.> Cit. Confusio #
# Author, Created: TP, May 2018 (version 18.1.0)
# Name: "Rabdomyzer"
# Description: "Un programmino per trovare le varianti"
# Updates: 18.1.1 (July 2018), 18.2.0 (August 2018), 21.1.0 (January 2021)
my $version="v21.1.0";

use warnings;
use strict;
use Data::Dumper qw(Dumper);
use List::MoreUtils qw(any uniq all);
use Getopt::Long qw(GetOptions);
use Path::Tiny;
use Pod::Usage qw (pod2usage);

                                            #Declaring command-line option arguments
my($vcf_file);
my($maf_threshold)=0.01;
my($model_file);
my($output_dir);
my($candidate);
my($gene_file);
my($help);

GetOptions( 'i=s' => \$vcf_file,
            'maf=f' => \$maf_threshold,
            'm=s' => \$model_file,
            'o=s' => \$output_dir,
            'gene=s' => \$candidate,
            'db=s' => \$gene_file,
            'h=s' => \$help
) or die "Invalid arguments! Enter $0 -h for more information";

die "Missing input VCF file specification" unless $vcf_file;
die "Missing output directory specification" unless $output_dir;
die "Missing model file specification" unless $model_file;

pod2usage(1) if $help;

=head1 NAME

Rabdomyzer $version

=head1 SYNOPSYS

$0 -i <file.vcf> -f 0.01 -m <file.model> -o <output.dir>

=head1 OPTIONS

-db    path to the gene annotation file
-gene  path to the candidate gene file
-f     minor allele frequency (default:0.01)

=head1 DESCRIPTION

The input file MUST be a VCF (Variant Call Format) file (for format specification, see: https://gatkforums.broadinstitute.org/gatk/discussion/11005/vcf-variant-call-format)
The input VCF file is annotated in the following way using Ensembl Variant Effect Predictor (VEP) v101:
vep --input_file <input.vcf> --output_file <output.VEP.vcf> --stats_file <output.VEP.stats.html> --warning_file <output.VEP.warning> --fork <N> --offline --force_overwrite --vcf --vcf_INFO_field ANN --fasta /opt/db/hg19/ucsc.hg19.fasta --species homo_sapiens --assembly GRCh37  --format vcf --everything --plugin SpliceRegion

Within the VCF file, the INFO field annotation section "ANN" MUST populated as described by the ANN.INFO.txt file

=cut

                                            #Declaring $ variables
my($i) = 0; my($j) = 0; my($k) = 0; my($x) = 0; my($y) = 0; my($z) = 0; my($ac) = 0;
my($candidate_prefix) = 0; my($most_sev) = 0;
my($ann_string) = 0; my($ann_csq_field) = 0; my($ann_maxaf_field) = 0; my($ann_canon_field) = 0; my($ann_biotype_field) = 0; my($ann_gene_field) = 0; my($ann_ccds_field) = 0;
my($case_geno) = 0; my($ctrl_geno) = 0;
my($key_gene) = 0; my($key_csq) = 0; my($key_onto) = 0; my($info_ann) = 0;
my($father) = 0; my($mother) = 0; my($parents) = 0; my($fath_col) = 0; my($moth_col) = 0; my($n_case) = 0; my($n_ctrl) = 0; my($n_unkn) = 0; my($n_geno) = 0; my($proband) = 0;
my($model) = 0;
my($n_var_fields) = 0;
my($geno_ctrl) = 0; my($geno_case) = 0; my($parent_field) = 0;
my($gnomad_var) = 0; my($mim_col) = 0; my($db_candidate) = 0; my($gene_name) = 0; my($name_candidate) = 0; my($parent_dir);

                                            #Declaring @ variables
my(@line) = (); my(@info) = (); my(@ann_var) = (); my(@ann) = (); my(@coord) = ();
my(@geno_var) = (); my(@iso) = (); my(@ctrl_ids) = (); my(@case_ids) = (); my(@unkn_ids) = ();
my(@ctrl_col) = (); my(@case_col) = (); my(@unkn_col) = (); my(@gene) = (); my(@format) = ();
my(@uniq_gene) = (); my(@uniq_inds_x_gene) = (); my(@ind_ids) = (); my(@fath_var_index) = (); my(@moth_var_index) = ();
my(@multi_csq) = (); my(@gene_header) = ();

                                            #Declaring/creating % variables
                                            
my(%iso_var) = ();
my(%gene_var) = ();
my(%inds_x_gene) = ();
my(%var_count) = ();
my(%gene_count) = ();
my(%gene_rep) = ();
my(%all_count) = ();
my(%db_gene) = ();

open ONTO, "$onto_file" or die $!;

while (<ONTO>)
 {
 next if ($_  =~ /^\#/);
 chomp;
 s/\r//;
 my $onto{$_}++;
 }
 
 close ONTO;
 
                                             # Checking/creating output folder

print "Checking output folder...\n\n";

if(-d "$output_dir")
 {
 print "Output directory already exists. No need to create it\n\n";
 }
 else
 {
 mkdir ("$output_dir");
 print "Created output directory\n";
 }

                                            # Saying what minor allele frequency is used for filtering

print "Using threshold for minor allele frequency: $maf_threshold\n\n";

                                            # Preparing gene database

if ($candidate)
   {
   $name_candidate = path($candidate)->basename;
   $parent_dir = path($output_dir)->parent; 
   print "Gene database will be placed in $parent_dir\n\n";
   $db_candidate = substr($name_candidate, 0, index($name_candidate, '.'));
   if (-f $parent_dir."/".$db_candidate.".db_gene.rabdomyzer")
      {
      print "Candidate gene database already exists. No need to prepare it\n\n";
      }
   else
      {
      print "Preparing gene database with candidate gene annotation... \n";
      system qq(perl rabdomyzer-match.pl -l $gene_file -q $candidate -o $parent_dir -R -H); #IMPORTANT: replace test.db_gene.rabdomyzer with abs path for dbNSFP3.5_gene.complete.rabdomyzer (e.g. /mnt/storage_slow/dbNSFP3.5a/dbNSFP3.5_gene.complete.rabdomyzer)
      print "Created file $db_candidate... Gene database prepared!\n\n"
      }
   open DBGENE, $parent_dir."/".$db_candidate.".db_gene.rabdomyzer" or die $!;
   }
else
   {
   print "Using default gene database dbNSFP3.5_gene.complete.rabdomyzer...\n\n";
   $db_candidate = "db";
   open DBGENE, "$gene_file" or die $!; #IMPORTANT: replace test.db_gene.rabdomyzer with abs path for dbNSFP3.5_gene.complete.rabdomyzer (e.g. /mnt/storage_slow/dbNSFP3.5a/dbNSFP3.5_gene.complete.rabdomyzer)
   }

while (<DBGENE>)
 {
 chomp;
 s/\r//;
 #s/^\s+//;
 #s/\s+$//;
 #next unless length;
 @line=split("\t", $_);
 if ($.==1)
    {
    @gene_header = @line;
    for $i ( 0..$#line )
       {
       if ($line[$i] eq "MIM_id")
          {
          $mim_col=$i;
          }
       }
    }
 else
    {
    if ($line[$mim_col] && $line[$mim_col] ne "None" && $line[$mim_col] ne "." && $line[$mim_col] =~ /^\d+$/)
       {
       $line[$mim_col]="https://www.omim.org/entry/".$line[$mim_col];
       }
    $gene_name = shift @line;
    $db_gene{$gene_name} = [ @line ];
    }
 }

shift @gene_header; 
                                                                                       
                                                                     # print Dumper($db_gene{MTOR})."\n\n\n";
                                                                     
                                            # Retrieving information from MODEL file

open MODEL, "$model_file" or die $!;

while (<MODEL>)
 {
 chomp;
 #s/\r//;
 s/#.*//;
 s/^\s+//;
 s/\s+$//;
 next unless length;
 @line=split("\t", $_);
 @case_ids=split(",", $line[0]);
 $n_case=scalar @case_ids;
 $proband=$case_ids[0];
 if ($line[3] ne "0")
    {
    @ctrl_ids=split(",", $line[3]);
    $n_ctrl=scalar @ctrl_ids;
    }
 if ($line[4] ne "0")
    {
    @unkn_ids=split(",", $line[4]);
    $n_unkn=scalar @unkn_ids;
    }    
 if ($line[1] ne "0" && $line[2] ne "0")
    {
    $model = "trio";
    $father = $line[1];
    $mother = $line[2];
    }
 elsif ($line[1] eq "0" && $line[2] eq "0")
    {
    if (scalar @case_ids==1 && $line[3] eq "0" && $line[4] eq "0")
       {
       $model = "single";
       }
    else
       {
       $model = "sharing";
       }
    }
 else
    {
    print "Failed to determine analysis model\n" and die $!;
    }
 
 print "Starting analysis of family $proband...\n\nModel in effect: $model\nCases are:$line[0] (Proband:$proband)\nControls are:$line[3]\nUnknowns are:$line[4]\n";
 if ($model eq "trio")
    {
    print "Proband's father and mother are:$father,$mother\n\n";
    }
 else
    {
    print "\n";
    }
 
 print "Creating tab-delimited file for heterozygous calls...\n";
 print "Creating tab-delimited file for homozygous calls...\n";

 open HET, "+>$output_dir/$proband.$model.$maf_threshold.$db_candidate.het_calls.txt" or die $!;
 open HOM, "+>$output_dir/$proband.$model.$maf_threshold.$db_candidate.hom_calls.txt" or die $!;
 if ($model eq "trio")
    {
    print "Creating tab-delimited file for compound heterozygous calls...\n\n";
    open COMHET, "+>$output_dir/$proband.$model.$maf_threshold.$db_candidate.comhet_calls.txt" or die $!;
    }
 else
    {
    print "\n";
    }
 }
 
print "Starting variant analysis...\n\n";
print "Performing consequence-based variant selection... \n\n";

open VCF, "$vcf_file" or die $!;

while (<VCF>)
 {
 chomp;s/\r//;
                                            # Retrieving information from header lines
 if ($_  =~ /^\#/)
    {
     if ($_  =~ /^\#\#INFO=<ID=ANN/)
       {
       $ann_string=substr($_, 99, -2);
                                            # Taking substring containing annotations separated by " | ". Substring positions may change with changing annotation. Must be checked every time annotation features change
       @ann=split("\\|", $ann_string);
       @ann=@ann[3,1,2,10,11,12,13,14,15,16,8,9,24,27,7,4,5,6,28,35,55,56,77,78,66,65,72,33,34,80,69,79,82,83,84,17,32,57,73,74,75,76,58,59,60,85,18,19,20,21,22,23,0,37];
       print "The following ".scalar @ann." variant-based annotation features are reported: Gene/Variant Number";
       for $i (0..$#ann)
          {
          print ",$ann[$i]";
          $ann_csq_field=$i if ($ann[$i] eq "Consequence");
          $ann_maxaf_field=$i if ($ann[$i] eq "MAX_AF");
          $ann_canon_field=$i if ($ann[$i] eq "CANONICAL");
          $ann_biotype_field=$i if ($ann[$i] eq "BIOTYPE");
          $ann_gene_field=$i if ($ann[$i] eq "SYMBOL");
          $ann_ccds_field=$i if ($ann[$i] eq "CCDS");
          }
       print "\n\n";      
       print "The following ".scalar @gene_header." gene-based annotation features are reported:";
       for $i (0..$#gene_header-1)
          {
          print "$gene_header[$i],";
          }
       print "$gene_header[$#gene_header]\n\n"; 
       }

    if ($_  =~ /^#CHROM/)                                                                # Exploring the sample information header line
      {
      @line=split("\t", $_);
      print HET "GENE/Var_N\t$line[0]\t$line[1]\t$line[3]\t$line[4]\t$line[6]\tGNOMAD_VAR\tAC";
      print HOM "GENE/Var_N\t$line[0]\t$line[1]\t$line[3]\t$line[4]\t$line[6]\tGNOMAD_VAR\tAC";
      if ($model eq "trio")
         {
         print COMHET "GENE/Var_N\t$line[0]\t$line[1]\t$line[3]\t$line[4]\t$line[6]\tGNOMAD_VAR\tAC";
         }
      
      for $i (0..$#ann)
         {
         print HET "\t$ann[$i]";
         print HOM "\t$ann[$i]";
         if ($model eq "trio")
            {
            print COMHET "\t$ann[$i]";
            }
         }   
      for $i (9..$#line)                                                                 # Looping on sample ids in VCF header
         {
         if (any {$_ eq $line[$i]} @case_ids)                                            # Checking what sample ids belong to "cases" and creating the case column index array
            {
            push (@case_col, $i);
            print HET "\t$line[$i]";
            print HOM "\t$line[$i]";
            if ($model eq "trio")
               {
               print COMHET "\t$line[$i]";
               }
            }
         }
      if ($model eq "trio")
         {
         for $i (9..$#line)
            {          
            if ($line[$i] eq $father)
               {
               $fath_col=$i;
               print HET "\t$line[$i]";
               print HOM "\t$line[$i]";
               print COMHET "\t$line[$i]";
               }
            }
         for $i (9..$#line)
            {          
            if ($line[$i] eq $mother)
               {
               $moth_col=$i;
               print HET "\t$line[$i]";
               print HOM "\t$line[$i]";
               print COMHET "\t$line[$i]";
               }
            }
         }
         for $i (9..$#line)
            {          
            if (any {$_ eq $line[$i]} @ctrl_ids)                                            # Checking what sample ids belong to "controls" and creating the control column index array
               {
               push (@ctrl_col, $i);
               print HET "\t$line[$i]";
               print HOM "\t$line[$i]";
               if ($model eq "trio")
                  {
                  print COMHET "\t$line[$i]";
                  }
               }
            }
         for $i (9..$#line)
            {          
            if (any {$_ eq $line[$i]} @unkn_ids)
               {
               push (@unkn_col, $i);                                                        # Checking what sample ids belong to "unknowns" and creating the case column index array
               print HET "\t$line[$i]";
               print HOM "\t$line[$i]";
               if ($model eq "trio")
                  {
                  print COMHET "\t$line[$i]";
                  }
               }   
            }
         for $i (0..$#gene_header)
            {
            print HET "\t$gene_header[$i]";
            print HOM "\t$gene_header[$i]";
            if ($model eq "trio")
               {
               print COMHET "\t$gene_header[$i]";
               }
            }                                                        # print "@ctrl_col\n"; # CHECKPOINT: array of control fields
         print HET "\n";
         print HOM "\n";
         if ($model eq "trio")
            {
            print COMHET "\n";
            }                                                        # print "@case_col\n"; # CHECKPOINT: array of case fields
         }                                                                                                                                                                 
    }                                                                                 

 
                                            # Looping through the variant sections
 else                                                                                    # Opening the "else" condition to parse the variant vcf section
 
    {
    @line=split("\t", $_);
    next if ($line[4] =~ /[\*\,]/);                                                      # Skipping multiallelic sites (ALT field with * or ,)              
    for $i (0..$#case_col)
       {
       @format=split(":", $line[$case_col[$i]]);
       push @geno_var, $format[0];
       }
    if ($model eq "trio")
       {
       @format=split(":", $line[$fath_col]);
       push @geno_var, $format[0];
       @format=split(":", $line[$moth_col]);
       push @geno_var, $format[0];
       }
    for $i (0..$#ctrl_col)
       {
       @format=split(":", $line[$ctrl_col[$i]]);         
       push @geno_var, $format[0];
       }
       for $i (0..$#unkn_col)
       {
       @format=split(":", $line[$unkn_col[$i]]);         
       push @geno_var, $format[0];
       }
       
    @info=split(";", $line[7]);                                                          # Creating array containing elements of the info field
                                                                                         # print "@info\n"; # CHECKPOINT: array of the info elements
    for $i (0..$#info)                                                                   # Looping through the info array elements
      {
      if ($info[$i] =~ /^AC\=/)
         {
         $ac = substr ($info[$i],3);
         }
      if ($info[$i] =~ /^ANN\=/)
         {
         $info_ann = substr ($info[$i],4);
         }
      else
         {
         $info_ann = $info[$i];
         }                                                                               # Removing the "ANN=" substring
      @ann_var=split(",", $info_ann);                                                    # Creating array containing annotation information for the different isoforms
      }
    
    for $i (0..$#ann_var)                                                                # Looping through the ann_var array
       {
       @iso=split("\\|", $ann_var[$i], -1);                                              # Creating array containing the annotation elements of the different isoforms in ann_var       
       next if !@iso;                                                                    # Checking whether @iso is empty (the array of the isoform is not populated at all)
       $iso[4]="http://gnomad.broadinstitute.org/gene/".$iso[4];                                                              
       @iso=@iso[3,1,2,10,11,12,13,14,15,16,8,9,24,27,7,4,5,6,28,35,55,56,77,78,66,65,72,33,34,80,69,79,82,83,84,17,32,57,73,74,75,76,58,59,60,85,18,19,20,21,22,23,0,37];
                                                                     # print "@iso\n"; # CHECKPOINT: array of the isoform elements
       if ( $iso[$ann_gene_field] && $iso[$ann_csq_field] && (!$iso[$ann_maxaf_field] || $iso[$ann_maxaf_field] eq "." || $iso[$ann_maxaf_field]<=$maf_threshold) && $iso[$ann_biotype_field] eq "protein_coding")                                                                                     # Verifying the following conditions: variant frequency <= frequency threshold, biotype = protein_coding
         {
         if ($iso[$ann_csq_field] =~ /\&/)
            {
            @multi_csq=split("&", $iso[$ann_csq_field]);
            foreach $key_onto (sort { $onto{$a} <=> $onto{$b} } keys %onto)
               {
               for $i (0..$#multi_csq)
                  {
                  if ($multi_csq[$i] eq $key_onto && $most_sev==0)
                     { 
                     splice (@iso, $ann_csq_field, 1, $multi_csq[$i]);
                     $most_sev++;
                     }
                  }     
               } 
            $most_sev=0;            
            }
        
         $iso_var{$iso[$ann_gene_field]}{$iso[$ann_csq_field]} = [ @iso ];
         if ($line[0] =~ /^chr[XYM\d]+$/)
            {
            $line[0] = substr $line[0], 3;
            $gnomad_var="http://gnomad.broadinstitute.org/variant/".$line[0]."-".$line[1]."-".$line[3]."-".$line[4];
            }
         else
            {
            $gnomad_var="http://gnomad.broadinstitute.org/variant/".$line[0]."-".$line[1]."-".$line[3]."-".$line[4];
            }
         @coord = ($line[0], $line[1], $line[3], $line[4], $line[6],$gnomad_var, $ac);
         unshift @{ $iso_var{$iso[$ann_gene_field]}{$iso[$ann_csq_field]} }, @coord;
         push @{ $iso_var{$iso[$ann_gene_field]}{$iso[$ann_csq_field]} }, @geno_var;
         }
      }

                                                                     # print Dumper(\%iso_var)."\n\n\n";

   foreach $key_gene (keys %iso_var)
      {
         foreach $key_onto (sort { $onto{$a} <=> $onto{$b} } keys %onto)
            {
            if (exists $iso_var{$key_gene}{$key_onto} && $most_sev==0)
               { 
               push @{$gene_var{$key_gene}}, $iso_var{$key_gene}{$key_onto};
               $most_sev++;
               }  
            } 
            $most_sev=0;
      }

   }

 $case_geno=0;
 $ctrl_geno=0;
 %iso_var=();
 @geno_var=();

 }                                                                                       # Closing loop on VCF file

print "Consequence-based variant selection done!\n\n";

$n_var_fields=scalar @ann + scalar @coord;
 
                                                                     # print Dumper(\%gene_var);
print "Starting model-based variant filtering... ";

 
 if ($model eq "trio")
    {
    $parents=2;
    foreach $key_gene (keys %gene_var)
       {
       for $i (0..$#{ $gene_var{$key_gene} } )
          {
          if ( ${$gene_var{$key_gene}}[$i][$n_var_fields+$n_case] =~ /^0[\|\/]0$/ &&  ${$gene_var{$key_gene}}[$i][$n_var_fields+$n_case+1] =~ /^0[\|\/]0$/ )
             {
             het( \@{$gene_var{$key_gene}[$i]}, \%db_gene );
             homdn( \@{$gene_var{$key_gene}[$i]}, \%db_gene );
             }
          if ( (${$gene_var{$key_gene}}[$i][$n_var_fields+$n_case] =~ /^0[\|\/]1$/ &&  ${$gene_var{$key_gene}}[$i][$n_var_fields+$n_case+1] =~ /^0[\|\/]0$/) || (${$gene_var{$key_gene}}[$i][$n_var_fields+$n_case] =~ /^0[\|\/]0$/ &&  ${$gene_var{$key_gene}}[$i][$n_var_fields+$n_case+1] =~ /^0[\|\/]1$/) || (${$gene_var{$key_gene}}[$i][$n_var_fields+$n_case] =~ /^0[\|\/]1$/ &&  ${$gene_var{$key_gene}}[$i][$n_var_fields+$n_case+1] =~ /^0[\|\/]1$/)  )
             {
             hom( \@{$gene_var{$key_gene}[$i]}, \%db_gene );
             }
          if ( ${$gene_var{$key_gene}}[$i][$n_var_fields+$n_case] =~ /^0[\|\/]1$/ && ${$gene_var{$key_gene}}[$i][$n_var_fields+$n_case+1] =~ /^0[\|\/]0$/ )
             {
             $parent_field=$i;
             push @fath_var_index, parentindex( @{$gene_var{$key_gene}[$i]} );
             }
          if ( ${$gene_var{$key_gene}}[$i][$n_var_fields+$n_case] =~ /^0[\|\/]0$/ &&  ${$gene_var{$key_gene}}[$i][$n_var_fields+$n_case+1] =~ /^0[\|\/]1$/ )
             {
             $parent_field=$i;
             push @moth_var_index, parentindex( @{$gene_var{$key_gene}[$i]} );
             }
          }
          for $i (0..$#fath_var_index)
             {
             for $k (0..$#moth_var_index)
                {
                for $j ($n_var_fields+$n_case+$parents..$n_var_fields+$n_case+$parents+$n_ctrl-1)
                   {
                   if ( (${$gene_var{$key_gene}}[$fath_var_index[$i]][$j] =~ /^0[\|\/]0$/ && ${$gene_var{$key_gene}}[$moth_var_index[$k]][$j] =~ /^0[\|\/]1$/) || (${$gene_var{$key_gene}}[$fath_var_index[$i]][$j] =~ /^0[\|\/]1$/ && ${$gene_var{$key_gene}}[$moth_var_index[$k]][$j] =~ /^0[\|\/]0$/) || (${$gene_var{$key_gene}}[$fath_var_index[$i]][$j] =~ /^0[\|\/]0$/ && ${$gene_var{$key_gene}}[$moth_var_index[$k]][$j] =~ /^0[\|\/]0$/) )
                      {
                      $ctrl_geno++;
                      }
                   }
                if ($ctrl_geno==$n_geno)
                   {
                   $x++;
                   print COMHET "$key_gene/$x";
                   foreach $j ( @{$gene_var{$key_gene}[$fath_var_index[$i]]} )
                      {
                      print COMHET "\t$j";
                      }
                   foreach $i (0..$#{ $db_gene{$key_gene} } )
                      {
                      print COMHET "\t$db_gene{$key_gene}[$i]"
                      }
                   print COMHET "\n";
                   print COMHET "$key_gene/$x";
                   foreach $j ( @{$gene_var{$key_gene}[$moth_var_index[$k]]} )
                      {
                      print COMHET "\t$j";
                      }
                   foreach $i (0..$#{ $db_gene{$key_gene} } )
                      {
                      print COMHET "\t$db_gene{$key_gene}[$i]"
                      }                      
                   print COMHET "\n";
                   }
             $ctrl_geno=0;
             }
          } 
        @moth_var_index = ();
        @fath_var_index = ();  
        $x=0;
       }
       
    }
    
 if ($model eq "sharing")
    {
    $parents=0;
    foreach $key_gene (keys %gene_var)
       {
       for $i (0..$#{ $gene_var{$key_gene} } )
          {
          het ( \@{$gene_var{$key_gene}[$i]}, \%db_gene );
          hom ( \@{$gene_var{$key_gene}[$i]}, \%db_gene );
          }
       }
    }
   
 if ($model eq "single")
    {
    foreach $key_gene (keys %gene_var)
       {
       for $i (0..$#{ $gene_var{$key_gene} } )
          {
          if (${$gene_var{$key_gene}}[$i][$n_var_fields] =~ /^0[\|\/]1$/)
             {
             $k++;
             print HET "$key_gene/$k";
             foreach $j (@{$gene_var{$key_gene}[$i]})
                {
                print HET "\t$j";
                }
             foreach $y (0..$#{ $db_gene{$key_gene} } )
                {
                if ($db_gene{$key_gene}[$y])
                   {
                   print HET "\t$db_gene{$key_gene}[$y]";
                   }
                else
                   {
                   print HET "\tNone";
                   }
                }
             print HET "\n";            
             }
          if (${$gene_var{$key_gene}}[$i][$n_var_fields] =~ /^1[\|\/]1$/)
             {
             $x++;
             print HOM "$key_gene/$x";
             foreach $j (@{$gene_var{$key_gene}[$i]})
                {
                print HOM "\t$j";
                }
             foreach $y (0..$#{ $db_gene{$key_gene} } )
                {
                if ($db_gene{$key_gene}[$y])
                   {
                   print HOM "\t$db_gene{$key_gene}[$y]";
                   }
                else
                   {
                   print HOM "\tNone";
                   }
                }
             print HOM "\n";            
                }   
          }
          $k=0;
          $x=0;
       }   
    }           

 close VCF;
 close DBGENE;
 close MODEL;
 close HET;
 close HOM;
 if ($model eq "trio")
 {
 close COMHET;
 }
 
print "Done!\n\n";

print "Writing xlsx file... ";
 
if ($model eq "trio")
   {
   system qq(python rabdomyzer-excel.py --het $output_dir/$proband.$model.$maf_threshold.$db_candidate.het_calls.txt --hom $output_dir/$proband.$model.$maf_threshold.$db_candidate.hom_calls.txt --comhet $output_dir/$proband.$model.$maf_threshold.$db_candidate.comhet_calls.txt --output $output_dir/$proband.$model.$maf_threshold.$db_candidate.xlsx);
   }
else
   {
   system qq(python rabdomyzer-excel.py --het $output_dir/$proband.$model.$maf_threshold.$db_candidate.het_calls.txt --hom $output_dir/$proband.$model.$maf_threshold.$db_candidate.hom_calls.txt --output $output_dir/$proband.$model.$maf_threshold.$db_candidate.xlsx);
   }
   
print "Done!\n\nAnanlysis of family $proband finished!\n";

                                            # Declaring subroutines

sub het                                     # Subroutine for de novo heterozygous (trio) and heterozygous (sharing) calls
   {
   my ($array_var, $hash_db) = @_;
   my @array_var = @{ $array_var };
   my %hash_db = %{ $hash_db };
   my $int_key_gene = $array_var[7];
   for $j ($n_var_fields+$n_case+$parents..$n_var_fields+$n_case+$parents+$n_ctrl-1)
      {                                                              #print "HET ctrl:$j\n";
      if ( $array_var[$j] && $array_var[$j] =~ /^0[\|\/]0$/)
         {
         $geno_ctrl++;
         }
      }
   for $j ($n_var_fields..$n_var_fields-1+$n_case)
      {                                                              #print "HET case:$j\n";
      if ( $array_var[$j] =~ /^0[\|\/]1$/ )
         {
         $geno_case++;
         }
      }                                                              #print "HET cases: genotypes $geno_case and number $n_case, controls: genotype $geno_ctrl and number $n_ctrl\n";
   if ($geno_case==$n_case && $geno_ctrl==$n_ctrl)
      {
      $x++;
      print HET "$int_key_gene/$x";
                                                                     #print "HET genotypes and numbers are OK\n ";
      for $j ( @array_var )
         {
         print  HET "\t$j";
         }
      foreach $k (0..$#{ $hash_db{$int_key_gene} } )
         {
         if ($hash_db{$int_key_gene}[$k])
            {
            print HET "\t$hash_db{$int_key_gene}[$k]";
            }
         else
            {
            print HET "\tNA";
            }
         }
      print HET "\n";
      } 
   $geno_ctrl=0;
   $geno_case=0;
   $x=0;   
   }
   
sub homdn                                   # Subroutine for de novo homozygous(trio) calls
   {
   my ($array_var, $hash_db) = @_;
   my @array_var = @{ $array_var };
   my $int_key_gene = $array_var[7];
   my %hash_db = %{ $hash_db };
   for $j ($n_var_fields+$n_case+$parents..$n_var_fields+$n_case+$parents+$n_ctrl-1)
      {
      if ( $array_var[$j] && $array_var[$j] =~ /^0[\|\/]0$/)
         {
         $geno_ctrl++;
         }
      }
   for $j ($n_var_fields..$n_var_fields-1+$n_case)
      {
      if ( $array_var[$j] && $array_var[$j] =~ /^1[\|\/]1$/ )
         {
         $geno_case++;
         }
      }
   if ($geno_case==$n_case && $geno_ctrl==$n_ctrl)
      {
      $x++;
      print HOM "$int_key_gene/$x";
      for $j ( @array_var )
         {
         print HOM "\t$j";
         }
      foreach $k (0..$#{ $hash_db{$int_key_gene} } )
         {
         if ($hash_db{$int_key_gene}[$k])
            {
            print HOM "\t$hash_db{$int_key_gene}[$k]";
            }
         else
            {
            print HOM "\tNA";
            }
         }
      print HOM "\n";
      }
   $geno_ctrl=0;
   $geno_case=0;
   $x=0;
   }

sub hom                                     #Subroutine for homozygous calls (trio and sharing)
   {
   my ($array_var, $hash_db) = @_;
   my @array_var = @{ $array_var };
   my $int_key_gene = $array_var[7];
   my %hash_db = %{ $hash_db };
   for $j ($n_var_fields+$n_case+$parents..$n_var_fields+$n_case+$parents+$n_ctrl-1)
      {                                                              #print "HOM ctrl:$j\n";
      if ( $array_var[$j] && ( $array_var[$j] =~ /^0[\|\/]0$/ || $array_var[$j] =~ /^0[\|\/]1$/ ) )
         {
         $geno_ctrl++;
         }
      }
   for $j ($n_var_fields..$n_var_fields-1+$n_case)
      {                                                              #print "HOM case:$j\n";
      if ( $array_var[$j] =~ /^1[\|\/]1$/ )
         {
         $geno_case++;
         }
      }
                                                                     #print "HOM cases: genotypes $geno_case and number $n_case, controls: genotype $geno_ctrl and number $n_ctrl\n";
   if ($geno_case==$n_case && $geno_ctrl==$n_ctrl)
      {
      $x++;
      print HOM "$int_key_gene/$x";
                                                                     #print "HOM genotypes and numbers are OK\n ";
      for $j ( @array_var )
         {
         print  HOM "\t$j";
         }               
      foreach $k (0..$#{ $hash_db{$int_key_gene} } )
         {
         if ($hash_db{$int_key_gene}[$k])
            {
            print HOM "\t$hash_db{$int_key_gene}[$k]";
            }
         else
            {
            print HOM "\tNA";
            }
         }
      print HOM "\n";
      } 
   $geno_ctrl=0;
   $geno_case=0;
   $x=0;
   } 

sub parentindex                             #Subroutine for indexing variants of parents for compound heterozygous calls (trio)
   {
   my @var_index;
   for $j ($n_var_fields..$n_var_fields-1+$n_case)
      {
      if ( $_[$j] =~ /^0[\|\/]1$/ )
         {
         $geno_case++;
         }
      }
      if ($geno_case==$n_case)
         {
         push @var_index, $parent_field;
         }
      $geno_case=0;
      return @var_index; 
   }
#   ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
   
