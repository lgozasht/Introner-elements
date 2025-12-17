#!/usr/bin/perl
$window = 0;

use English;
#  Note: for very large introns, this program causes large amounts of intron sequence to be stored
#  in memory, which may lead to memory problems.  The current solution to this is to 

%translate = ('TCA','S','TCC','S','TCG','S','TCT','S','TTC','F','TTT','F','TTA','L','TTG','L','TAC','Y','TAT','Y','TAA','X','TAG','X','TGC','C','TGT','C','TGA','X','TGG','W','CTA','L','CTC','L','CTG','L','CTT','L','CCA','P','CCC','P','CCG','P','CCT','P','CAC','H','CAT','H','CAA','Q','CAG','Q','CGA','R','CGC','R','CGG','R','CGT','R','ATA','I','ATC','I','ATT','I','ATG','M','ACA','T','ACC','T','ACG','T','ACT','T','AAC','N','AAT','N','AAA','K','AAG','K','AGC','S','AGT','S','AGA','R','AGG','R','GTA','V','GTC','V','GTG','V','GTT','V','GCA','A','GCC','A','GCG','A','GCT','A','GAC','D','GAT','D','GAA','E','GAG','E','GGA','G','GGC','G','GGG','G','GGT','G');

($#ARGV >= 1) || die "usage: HOMOLOGOUS_INTRON.PL blah.exons-introns blahblah.exons-introns.\n-window=length can be included anywhere on the commandline.\n\n";

$run = int (rand () * 1000000);

$fullstartover = $startover = 1;
@copy = @ARGV;

for ($i = $#ARGV; $i >= 0; $i--) {
    if ($ARGV[$i] eq "-nostartover") {
	$startover = 0;
	$fullstartover = 0;
	splice (@ARGV,$i,1);
    }
    elsif ($ARGV[$i] eq "-partstartover") {
	$startover = 0;
	splice (@ARGV,$i,1);
    }
    elsif ($ARGV[$i] eq "-firstversus") {
	$firstversus = 1;
	splice (@ARGV,$i,1);
    }
    elsif ($ARGV[$i] =~ /-orthologsfile=/) {
        $orth_group_file_name = $POSTMATCH;
	$startover = 0;
	$fullstartover = 0;	
	$use_provided_ortholog_list = 1;
        splice (@ARGV,$i,1);
    }
    elsif ($ARGV[$i] =~ /-all/) {
	$all = 1;
	splice (@ARGV, $i, 1);
    }
    elsif ($ARGV[$i] =~ /-filter/) {
	$filter = 1;
	unless ($window) { $window = 10 }
	unless ($cutoff) { $cutoff = 0.3 }
	splice (@ARGV, $i, 1);
    }
    elsif ($ARGV[$i] =~ /-window=/) {
        $window = $POSTMATCH;
        splice (@ARGV,$i,1);
    }
    elsif ($ARGV[$i] =~ /-cutoff=/) {
        $cutoff = $POSTMATCH;
        splice (@ARGV,$i,1);
    }
    elsif ($ARGV[$i] =~ /-allowgaps/) {
	$gaps = 1;
        splice (@ARGV,$i,1);
    }
}


if ($use_provided_ortholog_list) {
    open (IN, $orth_group_file_name);
    while (<IN>) {
	foreach $g (/\S+/g) {
	    $needed_genes{$g}++;
	}
    }
}



@sp = @ARGV;

#  Get data, make exint.pro files
foreach $sp (@sp) {
    ($root{$sp}) = $sp =~ /(?:.+\/|)(\S*?)(\.|$)/;
    open (IN, $sp);
    unless ($use_provided_ortholog_list) {
	open (OUT, ">$root{$sp}.exint.pro");
	print STDERR "Making $root{$sp}.exint.pro\n";
    }
    
    $/ = ">";
    <IN>;
    while (<IN>) {
	chomp;
	($name,$sequence) = /(\S+).*?\n(.+)/s;
	(!$use_provided_ortholog_list) || $needed_genes{$name} || next;
	($sequence =~ /\S/) || next;

	if ($sequence =~ /[\/\\]/) {
	    $sequence =~ s/(\S+)\s+(\S+)\s+(\S+)/$1$2$3/g;
	    @{$codons{$name}} = $sequence =~ /\S+/g;
	}
	else {
	    $sequence = join ("", $sequence =~ /\S+/g);
	    @{$codons{$name}} = $sequence =~ /[a-z]*[A-Z][a-z]*[A-Z][a-z]*[A-Z]/g;
	}



	$sequence =~ s/[[\/\\].+?[\/\\]//g;
	$sequence =~ s/ //g;
	while ($sequence =~ s/[^A-Z]+//) {
	    $intron_positions{$name} .= sprintf "%d.%d ", $-[0]/3+1,$-[0]%3;
	}
	$sequence =~ s/[A-Z]{3}/$translate{$&}/g;
	$pro{$name} = $sequence;
	$use_provided_ortholog_list && next;
	($pro{$name} =~ /\S/) && 
	    print OUT ">$name $intron_positions{$name}\n$pro{$name}\n";
	
    }
}


unless ($use_provided_ortholog_list) {
print STDERR "Out of the exint.pro section...\n";
#  Do all-against-all blast searches in both directions
foreach $sp1 (@sp) {
    system "formatdb -i $root{$sp1}.exint.pro";
    foreach $sp2 (@sp) {
	($sp1 eq $sp2) && next;
	print STDERR "In the blast section...\n";
	(open (IN, "$root{$sp2}.v.$root{$sp1}") && !$fullstartover) && next;
	print STDERR "blasting $root{$sp2}.exint.pro against $root{$sp1}.exint.pro\n";
	system "blastall.pl -i $root{$sp2}.exint.pro -d $root{$sp1}.exint.pro -p blastp -e0.0000000001 -m8 > $root{$sp2}.v.$root{$sp1}";
    }
}

#  Define ortholog pairs by best reciprocal hits

foreach $i (0..$#sp) {
    foreach $j ($i+1..$#sp) {
	print STDERR "Defining orthologs for $sp[$i] and $sp[$j].\n";
	push (@file_list, $outfile = join (".v.", sort ($root{$sp[$i]},$root{$sp[$j]})) . ".orthologs");
	unless ($use_provided_ortholog_list) {
	    system "perl ./ORTHOLOG_PAIRS.PL $root{$sp[$i]}.v.$root{$sp[$j]} $root{$sp[$j]}.v.$root{$sp[$i]} > $outfile";
	}
	push (@file_list, (join (".v.", sort ($root{$sp[$i]},$root{$sp[$j]})) . ".orthologs"));    
    }
}



close OUT;
}
($orth_group_file_name =~ /\S/) || 
    ($orth_group_file_name = join (".v.", sort values %root) . ".orthologs");


(@sp>2) &&
    ($startover || !open (IN, $orth_group_file_name)) &&
    system "ORTHOLOG_GROUPS.PL @file_list > $orth_group_file_name";



if ($all) {
    $orth_group_file_name = "HI.ListOfAllGenes";
    open (OUT, ">$orth_group_file_name");
    foreach $g (keys %pro) {
	print OUT "$g\t";
    }
    print OUT "\n";
    close OUT;
}
    


#  Align and report ortholog sets
$/ = "\n";
open (ORTH, "$orth_group_file_name");
while (<ORTH>) {
    $set_numb++;
    (@genes) = /\S+/g;
    open (OUT, ">clustin$run");
    foreach $gene (@genes) {
	print OUT ">$gene\n$pro{$gene}\n";
    }
    system "nohup clustalw clustin$run -output=gde -case=upper -outorder=input -quiet";
    open (IN, "clustin$run.gde");
    {
	local $/ = "%";
	$j = <IN>;
	$g = 0;
	while (<IN>) {
	    chomp;
	    $pos = 0;
	    ($seq) = /\n(.+)/s;
	    $name = $genes[$g];
	    $seq =~ s/\s//g;

	    @{$clust_aa{$name}} = $seq =~ /./g;
	    $seq =~ s/-/--->/g;
	    $seq =~ s/[A-Z]/$codons{$name}[$pos++].">"/eg;

	    @{$clust_nt{$name}} = $seq =~ /[^>]+/g;
	    $g++;
	}
    }
    print ">$set_numb:" . join (",", @genes);
    print "\nNumb:" . $#{$clust_aa{$genes[0]}} . "\n";


    %id = ();
    %conserved = ();
    unless ($all) {
	foreach $g1 (0..$#genes) {
	    foreach $g2 ($g1+1..$#genes) {
		for $i (0..$#{$clust_aa{$genes[0]}}) {
		    $t = $clust_aa{$genes[$g1]}[$i].$clust_aa{$genes[$g2]}[$i]; 
		    if ($t eq "--") {$symbol = "G"}
		    elsif ($t =~ /^-/) {$symbol = "F"}
		    elsif ($t =~ /-/) {$symbol = "H"}
		    elsif ($t =~ /(.)\1/) {$symbol = "1"}
		    else {$symbol = 0}		    
		    $id{$genes[$g1]}{$genes[$g2]}[$i] = 
			$id{$genes[$g2]}{$genes[$g1]}[$i] = $symbol;
		    
		}
	    }
	}
	
    }
    for $i (0..$#{$clust_aa{$genes[0]}}) {
	foreach $gene (@genes) {
	    print "$clust_aa{$gene}[$i]\t";
	}
	foreach $gene (@genes) {
	    print "$clust_nt{$gene}[$i]\t";
	}
	if (!$all && $window) {
	    $no = 0;
	    foreach $g1 (0..$#genes) {
		$firstversus && $g1 && last;
		$no && last;
		foreach $g2 ($g1+1..$#genes) {
		    $t = join ("", reverse @{$id{$genes[$g1]}{$genes[$g2]}}[0..$i-1]);
		    $t =~ /(.*?[01].*?){0,$window}[FGH]*/;
		    $left = reverse $&;
		    $pos = $id{$genes[$g1]}{$genes[$g2]}[$i];
		    $t = join ("", @{$id{$genes[$g1]}{$genes[$g2]}}[$i+1..$#{$clust_aa{$genes[0]}}]);
		    $t =~ /(.*?[01]){0,$window}[FGH]*/;
		    $right = $&;

		    if ($filter) {
			unless ($gaps || !(($left.$right) =~ /[FGH]/)) { $no++; last; }
			unless (($left =~ tr/01/01/) && ($right =~ tr/01/01/)) { $no++; last; }
			unless ($cutoff <= (($left =~ tr/1/1/) / ($left =~ tr/01/01/))) { $no++; last; }
			unless ($cutoff <= (($right =~ tr/1/1/) / ($right =~ tr/01/01/))) { $no++; last; }
		    }
		    else {
			print "$left($pos)$right:";
		    }


		}
	    }
	    if ($filter) {
		print (("Y","N")[$no > 0]);
	    }
	    

	}
	print "\n";

    }

}



	
