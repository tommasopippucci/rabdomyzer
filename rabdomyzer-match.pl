#!/usr/bin/env perl
                                                                                                   # < When family, no gene. When gene, no family.> Cit. Confusio #
use warnings;
use strict;
use Getopt::Long qw(GetOptions);
use Path::Tiny;

# Author: TP Nov 2018, modified from previous versions Dec 2014, Oct 2017
# Name: match_features.pl
# Aim: to match elements in LIST FILE and QUERY FILE according to content of fields LIST KEY in LIST FILE and QUERY KEY in QUERY FILE. Reports all QUERY_LIST and QUERY_FILE rows which did/did not match. 
   #Each row will be reported with LIST/QUERY KEY as first column.
# Usage: perl match_features.pl -l <LIST_FILE> -q <QUERY_FILE> -p <MATCH_PREFIX> -o <OUTPUT_DIR>
# Optional arguments:
# -lcol LIST_KEY (Default: 1)
# -qcol QUERY_KEY (Default: 1)
# -ncol Number of columns forming the matching element (Default: 1)
# -R Reciprocal match (Default: NO)
# -H If header is present (Default: NO). IMPORTANT: Header must be either present or absent in both LIST and QUERY files

                                            #Declaring option arguments

my($listin);
my($querin);
my($list_elem) = 0;
my($query_elem) = 0;
my($elem_length) = 1;
my($header);
my($reciprocal);
my($output_dir);

GetOptions( 'l=s' => \$listin,
            'q=s' => \$querin,
            'lcol=i' => \$list_elem-1,
            'qcol=i' => \$query_elem-1,
            'ncol=i' => \$elem_length,
            'H' => \$header,
            'R' => \$reciprocal,
            'o=s' => \$output_dir,         
) or die "Invalid arguments!";

die "Missing list file specification" unless $listin;
die "Missing query file specification " unless $querin;
die "Missing output directory specification" unless $output_dir;

if(-d "$output_dir")
 {
 #print "Output directory already exists. No need to create it\n\n";
 }
 else
 {
 mkdir ("$output_dir");
 #print "Created output directory\n";
 }
                                            #Declaring $,@,% variables
my($elem_pos) = 0;
my(@complete_header) = ();
my(%list) = ();
my(@line) = ();
my($query_match) = 0;
my(@list_header) = ();
my(@query_header) = ();
my($elements) = 0;
my($list_key) = 0;
my($match) = 0;
my(@match) = ();
my($match_count) = 0;
my($i) = 0;
my(@decomp) = ();
my $name_querin = path($querin)->basename;
my $matchprefix = substr($name_querin, 0, index($name_querin, '.'));

                                            #Printing script description to STDOUT
open OUT, "+>$output_dir/$matchprefix.db_gene.rabdomyzer" or die $!;

open LIST, "$listin" or die $!;
                                            #Loading LIST header, keys and elements
while (<LIST>)
 {
 chomp;
 s/\r//;
 @line=split('\t', $_, -1);
 $elem_pos=$list_elem;
 match_elem_arg(@line);
 if (defined $header and $. == 1)
    {
    @list_header = @line;
    next;
    }
 $list {$match} = [ @line ];
 $elements++;
 }

print "   Elements in $listin: $elements...\n";
$elements=0;
                                            #Performing LIST-QUERY match
open QUERY, "$querin" or die $!;

while (<QUERY>)
 {
 next if $_ eq "\n"; 
 chomp;
 s/\r//;
 @line=split('\t', $_, -1);
 $elem_pos=$query_elem;
 match_elem_arg(@line);
                                            #Composing output header
 if (defined $header && $. == 1)
    {
    @query_header = @line;
    @complete_header = (@query_header, @list_header);
    print_match($match);
    for ( $i=0 ; $i<=$#complete_header ; $i++)
       {
       print OUT "\t$complete_header[$i]";
       }
    print OUT "\n";
    next;
    }
$elements++;
                                            #Instructions for when a match is found
 foreach $list_key (keys %list)
    { 
    if ( $match eq $list_key) 
       {
       print_match($match);
       for $i (0 .. $#line)
          {
          print OUT "\t$line[$i]";
          }
       for $i ( 0 .. $#{ $list{$match} } )
          {
          print OUT "\t$list{$match}[$i]";
          }
       print OUT "\n";   
       $query_match++;
       $match_count++;
       if (defined $reciprocal)
          {
          delete $list{$list_key};
          }
       }
   }      
                                            #Instructions for when there is no matching key in LIST
 if ($query_match==0)
    {
    print_match($match);
    for ($i=0 ; $i<=$#line ; $i++)
       {
       print OUT "\t$line[$i]";
       }
    for ($i=0 ; $i<=$#list_header ; $i++)
       {   
       print OUT "\tNA";
       }
    print OUT "\n";
   }
 $query_match=0;
 
 }

print "   Elements in $querin: $elements...\n";
print "   Matching elements: $match_count...\n";
if (defined $reciprocal)
   {
   print "   Elements unique to $querin: ";
   print $.-$match_count-1;
   print "...\n";
   }

                                            #Instructions for when there is no matching key in QUERY

$match_count=0;

if (defined $reciprocal)
   {
   foreach $list_key (sort {lc $a cmp lc $b } keys %list)
      {
      $match_count++;
      print_match($list_key);
      for ( $i=0 ; $i<=$#query_header ; $i++)
         {   
         print OUT "\tNA";
         }
      for $i ( 0 .. $#{ $list{$list_key} } )
         {
         print OUT "\t$list{$list_key}[$i]";
         }
   print OUT "\n";
      }
print "   Elements unique to $listin: $match_count...\n";
   }


close LIST;
close QUERY;
close OUT;

                                                           #Declaring subroutines         

                                            #Declaring subroutine to create line-by-line match element and argument
               
sub match_elem_arg{
                  @match = splice @_, $elem_pos, $elem_length;
                  $match=$match[0];
                  if ($elem_length>1)
                     {
                     for $i (1..$#match)
                        {
                        $match=$match.":".$match[$i];
                        }
                     }
                  @line = @_;   
                  return @_;
                  return $match;                     
                  }

                                            #Declaring subroutine for printing simple/complex match elements made of 1/>1 original column fields
               
sub print_match{
                if ($elem_length==1)
                   {
                   print OUT $_[0];
                   }
                   else
                   {
                   @decomp = split(":", $_[0]);
                   for $i (0..$#decomp-1)
                      {
                      print OUT "$decomp[$i]\t";
                      }
                   print OUT $decomp[$#decomp];
                   }
               }

