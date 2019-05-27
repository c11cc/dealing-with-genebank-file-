#annotate the bla result to genes,work on this laptop
use strict;
use warnings;

print "enter genebank file:\t";
my $file=<>;

my %feature;
my $i=-1;
my $j=-1;
my $notes=0;
my $joi=0;
my @num;
my $mrn=0;

open IN, "$file" or die $!;
while(<IN>)
  {
    chomp;
    my $tmp=$_;    
    if($tmp=~/FEATURES/)
      {	
        if($j gt -1 and $i gt -1)
          {
            $num[$j]=$i;
          } 
        $j++;
        $i=-1;
        if($notes)
          {
            $notes=0;	
          } 
        if($joi)
          {
            $joi=0;	
          }                  
      }
    elsif($tmp=~/\/organism\=\"([^\t]+)\"/)#/organism="Yarrowia lipolytica"
      {
        $feature{organism}->[$j]=$1;      
        next;
      }      
    elsif($tmp=~/\/strain\=\"([^\t]+)\"/)#/strain="CLIB89(W29)"
      {
        $feature{strain}->[$j]=$1;   
        next;  
      }      
    elsif($tmp=~/\/mol_type\=\"([^\t]+)\"/)#/mol_type="genomic DNA"
      {
        $feature{mol_type}->[$j]=$1;  
        next;   
      } 
    elsif($tmp=~/\/organelle\=\"([^\t]+)\"/)#/organelle="mitochondrion"
      {
        $feature{organelle}->[$j]=$1; 
        next;    
      }
    elsif($tmp=~/\/chromosome\=\"([^\t]+)\"/)#/chromosome="1A"
      {
        $feature{chromosome}->[$j]=$1;     
        next;
      }
    elsif($tmp=~/VERSION[^\w]+([^\t]+)/)#VERSION     CP017553.1
      {
        $feature{chromosome_id}->[$j+1]=$1;     
        next;
      }      
     elsif($tmp=~/source\D+(\d\.\.\d+)$/)#source          1..2257857
      {
        $feature{source}->[$j]=$1;    
        next;
      }      
      
    elsif($tmp=~/gene[^\d]+(\d+\.\.\d+)/)#gene            3285..3920  /   gene            <13383..14354
      {
      	$i++;
       	$feature{start}->[$j][$i][0]=0;#whole region
      	$feature{start}->[$j][$i][1]=$1;#3285..3920 
      	if($notes)
          {
            $notes=0;	
          }  
        if($joi)
          {
            $joi=0;	
          }   	
      }
    elsif($tmp=~/gene[^\(]+\((\d+\.\.\d+)\)/)#gene            complement(1230..1648)
      {
      	$i++;
       	$feature{start}->[$j][$i][0]=0;#whole region
      	$feature{start}->[$j][$i][1]=$1;#1230..1648 
      	$feature{reverse}->[$j][$i]=1;
      	if($notes)
          {
            $notes=0;	
          } 
        if($joi)
          {
            $joi=0;	
          }    	
      }

    elsif($tmp=~/misc_feature[^\w]+(\d+\.\.\d+)/)#misc_feature    1635034..1635310
      {
      	$i++;
      	$feature{start}->[$j][$i][0]=0;# have no intron
      	$feature{start}->[$j][$i][1]=$1;#1230..1457,1523..1648
      	$feature{type}->[$j][$i]="other";#product type
      	if($notes)
          {
            $notes=0;	
          } 
        if($joi)
          {
            $joi=0;	
          }
      }   
    elsif($tmp=~/misc_feature[^\(]+\(([^\t]+)\)/)#misc_feature    complement(1859372..1859643)
      {
      	$i++;
      	$feature{start}->[$j][$i][0]=0;# have no intron
      	$feature{start}->[$j][$i][1]=$1;#1230..1457,1523..1648
      	$feature{type}->[$j][$i]="other";#product type
      	$feature{reverse}->[$j][$i]=1;
      	if($notes)
          {
            $notes=0;	
          } 
        if($joi)
          {
            $joi=0;	
          }
      }         
               
    elsif($tmp=~/CDS[^\(]+\(join\(([^\)]+)\)\)/)#CDS            complement(join(1230..1457,1523..1648)) 
      { 
      	$feature{start}->[$j][$i][0]=1;#have intron
      	$feature{start}->[$j][$i][2]=$1;#1230..1457,1523..1648
      	$feature{reverse}->[$j][$i]=1;      	
      	if($notes)
          {
            $notes=0;	
          }
        if($joi)
          {
            $joi=0;	
          } 
      }         
    elsif($tmp=~/CDS[^\(]+\(join\(([^\)]+)/)#CDS           complement(join(315569..315894,315955..315975,
      {
      	$feature{start}->[$j][$i][0]=1;#have intron
      	$feature{start}->[$j][$i][2]=$1;#1230..1457,1523..1648
      	$feature{reverse}->[$j][$i]=1;      	
      	if($notes)
          {
            $notes=0;	
          }
        $joi=1;
      }
    elsif($tmp=~/CDS\s+join\(([^\)]+)/)#CDS            join(315569..315894,315955..315975,
      {
      	$feature{start}->[$j][$i][0]=1;#have intron
      	$feature{start}->[$j][$i][2]=$1;#1230..1457,1523..1648      	
      	if($notes)
          {
            $notes=0;	
          }
        $joi=1;
      }           
    elsif($tmp=~/CDS[^\d]+(\d+\.\.\d+)/)#CDS            3285..3920#using ^\d but not \w may cause error for the complement instance
      {
      	$feature{start}->[$j][$i][0]=0;# have no intron
      	$feature{start}->[$j][$i][2]=$1;#3285..3920

      	if($notes)
          {
            $notes=0;	
          }
      }     
     elsif($tmp=~/tRNA[^\(]+\(([^\t]+)\)/)#tRNA            complement(231750..231834)
      {
      	if(!$i+1)
      	  {
      	    next;	
      	  }
      	$feature{start}->[$j][$i][0]=0;#have intron
      	$feature{start}->[$j][$i][2]=$1;#231750..231834
      	$feature{type}->[$j][$i]="RNA";#product type
      	$feature{reverse}->[$j][$i]=1;      	
      	if($notes)
          {
            $notes=0;	
          } 
        if($joi)
          {
            $joi=0;	
          }
      }   
    elsif($tmp=~/tRNA[^\w]+(\d+\.\.\d+)/)#tRNA            343716..343788
      {
      	$feature{start}->[$j][$i][0]=0;# have no intron
      	$feature{start}->[$j][$i][2]=$1;#1230..1457,1523..1648
      	$feature{type}->[$j][$i]="RNA";#product type
      	if($notes)
          {
            $notes=0;	
          } 
        if($joi)
          {
            $joi=0;	
          }
      }
      
     elsif($tmp=~/rRNA[^\(]+\(([^\t]+)\)/)#rRNA            complement(3136045..3136168)
      {
      	if(!$i+1)
      	  {
      	    next;	
      	  }
      	$feature{start}->[$j][$i][0]=0;#have intron
      	$feature{start}->[$j][$i][2]=$1;#231750..231834
      	$feature{type}->[$j][$i]="RNA";#product type
      	$feature{reverse}->[$j][$i]=1;      	
      	if($notes)
          {
            $notes=0;	
          } 
        if($joi)
          {
            $joi=0;	
          }
      }   
    elsif($tmp=~/rRNA[^\w]+(\d+\.\.\d+)/)#rRNA            3429182..3429300
      {
      	$feature{start}->[$j][$i][0]=0;# have no intron
      	$feature{start}->[$j][$i][2]=$1;#1230..1457,1523..1648
      	$feature{type}->[$j][$i]="RNA";#product type
      	if($notes)
          {
            $notes=0;	
          } 
        if($joi)
          {
            $joi=0;	
          }
      }
     elsif($tmp=~/ncRNA[^\(]+\(([^\t]+)\)/)#ncRNA           complement(1021230..1021309)
      {
      	if(!$i+1)
      	  {
      	    next;	
      	  }
      	$feature{start}->[$j][$i][0]=0;#have intron
      	$feature{start}->[$j][$i][2]=$1;#231750..231834
      	$feature{type}->[$j][$i]="RNA";#product type
      	$feature{reverse}->[$j][$i]=1;      	
      	if($notes)
          {
            $notes=0;	
          } 
        if($joi)
          {
            $joi=0;	
          }
      }   
    elsif($tmp=~/ncRNA[^\w]+(\d+\.\.\d+)/)#ncRNA           1370904..1371346
      {
      	$feature{start}->[$j][$i][0]=0;# have no intron
      	$feature{start}->[$j][$i][2]=$1;#1230..1457,1523..1648
      	$feature{type}->[$j][$i]="RNA";#product type
      	if($notes)
          {
            $notes=0;	
          } 
        if($joi)
          {
            $joi=0;	
          }         
      }      
      
    elsif($tmp=~/mRNA\s+/)#mRNA            complement(join(1230..1457,1523..1648))
      {
      	$feature{type}->[$j][$i]="protein";#product type 
        if($notes)
          {
            $notes=0;	
          }  
        if($joi)
          {
            $joi=0;	
          }    
        $mrn=1;      	    	
      }
     
    elsif($tmp=~/(\w+)\s+(\d+\.\.\d+)/)
      {
      	$i++;
      	$feature{start}->[$j][$i][0]=0;# have no intron
      	$feature{start}->[$j][$i][1]=$2;#
        $feature{type}->[$j][$i]=$1;
        $notes=0 if($notes);
        $joi=0 if($joi);
      }
           
    elsif($tmp=~/\/locus_tag\=\"([^\t]+)\"/)#/locus_tag="YALI1_A00014g"
      {
        $feature{name}->[$j][$i]=$1;
        if($notes)
          {
            $notes=0;	
          }  
        if($joi)
          {
            $joi=0;	
          }    
        $mrn=0;  
      }      
    elsif($tmp=~/\/product=\"([^\t]+)\"/)#/product="hypothetical protein"
      {
        $feature{product}->[$j][$i]=$1;
        if($notes)
          {
            $notes=0;	
          }
        if($joi)
          {
            $joi=0;	
          }
      }   
    elsif($tmp=~/\/protein_id\=\"([^\t]+)\"/)#/protein_id="AOW00072.1"
      {
        $feature{protein_id}->[$j][$i]=$1;
      }  
    elsif($tmp=~/\/standard_name\=\"([^\t]+)\"/)#/standard_name="ha2747"
      {
        $feature{protein_id}->[$j][$i]=$1;
      }     
    elsif($tmp=~/\/gene\=\"([^\t]+)\"/)#/gene="trnA(tgc)"
      {
        $feature{protein_id}->[$j][$i]=$1;
      }  
    elsif($tmp=~/pseudo/)#/gene="trnA(tgc)"
      {
        $feature{type}->[$j][$i].="pseudo";
      }
    elsif($tmp=~/\/note\=\"([^\"]+)/)# /note="Compare to YALI0A00110g, weakly similar to
      {
        $feature{note}->[$j][$i]=$1;
        $notes=1;
        $mrn=0;
        next;
      } 
    
    elsif($tmp=~/ORIGIN/)
      {
        $notes=0;	
      }
    elsif($notes and $tmp=~/\/codon_start/)
      {
        $notes=0;        
      }   
            
    elsif($notes eq 1 and $tmp=~/                     [^\t]+\s\w+/ and $tmp!~/(\/locus_tag)|(\/codon_start)|(\/product)|(\/protein_id)|(\/translation)|(\/gene)|(\/inference)/)
      { 
      	if($tmp=~/                     ([^\"]+)$/)
      	  {
            $feature{note}->[$j][$i].=" ".$1;     
          }   
        elsif($tmp=~/                     ([^\t]+)\"/)
      	  {
            $feature{note}->[$j][$i].=" ".$1;      
          }
      }     
                                                               # mRNA            complement(join(315569..315894,315955..315975,
    elsif($joi eq 1 and $tmp=~/[^\d]+([^\)]+)\)/)  #                     316049..316064))
      {
        $feature{start}->[$j][$i][2].=$1; 
        $joi=0;       
      }  
                                                                # mRNA            complement(join(315569..315894,315955..315975,
    elsif($joi eq 1 and $tmp=~/[^\d]+([^\)]+)$/)  #                     316049..316064,
      {
        $feature{start}->[$j][$i][2].=$1;        
      }                 
        
    if($_=~/complement/)
      {
      	$feature{reverse}->[$j][$i]=1;
      }          
  } 
close IN;
$num[$j]=$i;#the last one 
my $chro_num=$j;

open OUT, ">feature.re " or die $!;
print OUT "strain\torganelle\tgenome\tsource\tfunction_unit_type\tfunction_unit_id\tchromosome\tfunction_unit site_easy\tfunction_unit site_complex\tstrand\tfunction_unit_function\tgene_name\tfunction_unit_note\n";
foreach $j(0..$j)
  {
  	foreach $i(0..$num[$j])
  	  { 
        if($feature{organism}->[$j])
          {
          	print OUT "$feature{organism}->[$j]\t";
          }
        else
          {
            print OUT "na\t";
          }      
        if($feature{organelle}->[$j])
          {
          	print OUT "$feature{organelle}->[$j]\t";
          }
        else
          {
            print OUT "na\t";
          }  
        if($feature{chromosome_id}->[$j])
          {
          	print OUT "$feature{chromosome_id}->[$j]\t";
          }
        else
          {
            print OUT "na\t";	
          }
        if($feature{source}->[$j])
          {
          	print OUT "$feature{source}->[$j]\t";
          }
        else
          {
            print OUT "\t";
          } 
        if($feature{type}->[$j][$i])#product
          {
          	print OUT "$feature{type}->[$j][$i]\t";
          }
        else
          {
            print OUT "mit_protein\t";
          }                   
        if($feature{name}->[$j][$i])#product_id
          {
          	print OUT "$feature{name}->[$j][$i]\t";
          }
        else
          {
            print OUT "na\t";
          }           
        #if($feature{strain}->[$j])
        #  {
        #  	print OUT "$feature{strain}->[$j]\t";
        #  }
        #else
        #  {
        #    print OUT "na\t";
        #  }
        #if($feature{mol_type}->[$j])
        #  {
        #  	print OUT "$feature{mol_type}->[$j]\t";
        #  }
        #else
        #  {
        #    print OUT "na\t";
        #  }
        if($feature{chromosome}->[$j])#chromosome
          {
          	print OUT "$feature{chromosome}->[$j]\t";
          }
        else
          {
            print OUT "na\t";
          }
        if($feature{start}->[$j][$i][0])#site_easy/tsite complex
          {
          	print OUT "$feature{start}->[$j][$i][1]\t$feature{start}->[$j][$i][2]\t";
          }
        else
          {
            print OUT "$feature{start}->[$j][$i][1]\tna\t";
          } 
        if($feature{reverse}->[$j][$i]) #strand
          {
            print OUT "-\t";	
          }
        else
          {
            print OUT "+\t";	
          }
        
        if($feature{product}->[$j][$i])#product
          {
          	print OUT "$feature{product}->[$j][$i]\t";
          }
        else
          {
            print OUT "na\t";
          }    
        if($feature{protein_id}->[$j][$i])#product id
          {
          	print OUT "$feature{protein_id}->[$j][$i]\t";
          }
        else
          {
            print OUT "na\t";
          }      
        if($feature{note}->[$j][$i])#note
          {
          	print OUT "$feature{note}->[$j][$i]\t";
          }
        else
          {
            print OUT "na\t";
          }
        print OUT "\n";
      }
  }
close OUT;