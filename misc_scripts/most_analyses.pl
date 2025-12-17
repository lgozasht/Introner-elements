use English;
$maxlength = 100;

#  Section1 : Get relationships from species tree (first argument)
$treefile = $ARGV[0];
($treepointer,$taxapointer) = &maketree ($treefile,"$treefile.CladesAndDivergences");
%divergences = %$treepointer;
%alltaxa = %$taxapointer;
@roots = sort keys %alltaxa;


&checkforfiles(@roots);

#  Section2: Make exons-introns files, protein files, and diamond database files
#  Note, for each species, there must exist two files:
#  SpeciesName.gff (gff or gtf file, some flexibility of formatting)
#  SpeciesName.fna (genome fasta file)
foreach $root (@roots) {

    &makeexonsintrons ("$root.gff","$root.fna","$root.exons-introns"); # exons-introns file in $root.exons-introns
    &makeintronsflanks("$root.exons-introns","$root.intronsflanks"); # introns+flanking sequence file in $root.intronsflanks

    &makepro("$root.exons-introns","$root.pro"); # protein file in $root.pro
    &makedmnd("$root.pro"); # make diamond database file for $root.pro
}


# Section 3: Define orthologs by pairwise recipropcal blast, then find shared introns
for $r1 (@roots) {

    for $r2 (@roots) {
	last if $r1 eq $r2;
	&makeorthologs("$r1.pro","$r2.pro","$r1-$r2.orthologs"); # pairwise ortholog file in $r1-$r2.orthologs
	&makesharedintrons("$r1.exons-introns","$r2.exons-introns","$r1.pro","$r2.pro","$r1-$r2.orthologs","$r1-$r2.sharedintrons"); # shared introns file ("-" means unshared intron at well-aligned position) in $r1-$r2.sharedintrons
    }
}

#  Section 4: For each species, compile shared introns, then generate files
for $root (@roots) {

    &compilesharedintrons($root,\@roots,".sharedintrons","$treefile.CladesAndDivergences","$root.intronsharingcompiled"); # Compilation of intron sharing patterns in $r.intronsharingcompiled
    &makeintronsflankswithsharing("$root.intronsflanks","$root.intronsharingcompiled","$root.intronsflankssharing"); # File of introns+flanking sequence, with intron sharing positions indicated, in $root.intronsflankssharing
}


#  Section 5: Identify groups
#foreach $root (@roots) {
 #   &makeends("$root.intronsflankssharing","$root.ends");

#    &blastends("$root.ends","$root.ends.blastout");
#    &makeparalogs("$root.pro","$root.paralogs");
#    &filterblastends("$root.ends.blastout","$root.paralogs", "$root.ends.blastout.filtered");
#    &makematches ("$root.ends.blastout.filtered","$root.intronmatches");
#    %clusters = &makeclusters("$root.intronmatches");
#    %clusters = &orderbynumber(%clusters);
#    &printclusters ("$root.intronsflankssharing", "$root.intronflankssharing.IntronerGroups","$root.IntronersByGroup", %clusters);
#}

# Section 6: Build Malin matrix
&makeorthologsmultispecies(\@roots,"$treefile.orthologs");
$positions = &makesharedintronsmanypaml(\@roots,"$treefile.orthologs","$treefile.sharedintrons","$treefile.pamlin");
&makemalinmatrix(\@roots, "$treefile.pamlin", "$treefile.sharedintrons","$treefile.malin");

#  Section 7:
#  Phase and estimated phase, per group
#  5' dinucleotide, by group
#  5' dinucleotide, by phylogenetic distribution
#  evidence for cutter, by group
#  evidence for cutter, by phylogenetic distribution
#  evidence for cutter, by phylogenetic distribution, previously determined categories (for figure)

for $root (@roots) {
    $tetranucphasedistpointer = &maketetranuc("$root.exons-introns");
    &makeexpobsphase("$root.intronflankssharing.IntronerGroups",$tetranucphasedistpointer,"$root.ExpObsPhase");
    &donorss("$root.intronflankssharing.IntronerGroups","$root.DonorByPhylodist",5);
    &donorss("$root.intronflankssharing.IntronerGroups","$root.DonorByGroup",6)
    &cutter4("$root.intronflankssharing.IntronerGroups","$root.CutterByPhylodist",5);
    &cutter4("$root.intronflankssharing.IntronerGroups","$root.CutterByGroup",6);
}




sub donorss () {
    my ($infile,$outfile,$column) = @_;
    my %sstally = ();
    my $tot;
    my $head;
    local $/ = ">";
    open (IN, $infile);
    while (<IN>) {
	($head) = /(.+)/;
	$group = (split (/:/, $head))[$column-1];
	/\n[A-Z]+([a-z]{2})/;
	$sstally{$group}{$1}++;
    }
    open (OUT, ">$outfile");
    print OUT "Group\tGA\tGC\tGT\tTotal\n";
    for $group (sort keys %sstally) {
	$tot = $sstally{$group}{"ga"} + $sstally{$group}{"gc"} + $sstally{$group}{"gt"};
	next unless $tot;
	print OUT "$group\t";
	
	for $ss ("ga","gc","gt") {
	    printf OUT "%d\t", $sstally{$group}{$ss};
	}
	print OUT "$tot\n";
    }
    close OUT;
}

sub cutter4 () {
    my ($infile,$outfile,$column) = @_;
    my %cuttertally = ();
    my $tot;
    my $head;
    my $i;
    my $matches;
    local $/ = ">";
    open (IN, $infile);
    while (<IN>) {
	($head) = /(.+)/;
	$group = (split (/:/, $head))[$column-1];

	($donor,$fivesix,$downex) = /\n[A-Z]+([a-z]{4})([a-z]{2})[a-z]+([A-Z]{4})/;
	@donor = $donor =~ /./g;
	@downex = $downex =~ /./g;
	$matches = 0;
	for $i (0..3) {
	    $matches += $downex[$i] eq uc($donor[$i]);
	}
	$cuttertally{$group}{"$matches-$fivesix"}++;
    }

    open (OUT, ">$outfile");
    print OUT "Group\t4-ct\t4-cc\t3-ct\t3-cc\tTotal\n";
    for $group (sort keys %cuttertally) {
	print OUT "$group\t";
	$tot = sum (values %{$cuttertally{$group}});
	for $patt ("4-ct","4-cc","3-ct","3-cc") {
	    printf OUT "%d\t", $cuttertally{$group}{$patt};
	}
	print OUT "$tot\n";
    }
    close OUT;
}


sub orderbynumber () {
    my %clusters = @_;
    my %tally = ();
    my %convert = ();
    my $v;
    my $i;
    for $v (values %clusters) {
	$v =~ /\d+/;
	$tally{$&}++;
    }
    my @a = reverse sort {$tally{$a} <=> $tally{$b}} keys %tally;
    for $i (0..$#a) {
	$convert{$a[$i]} = $i+1;
    }
    for $i (keys %clusters) {
	$clusters{$i} =~ s/\d+/$convert{$&}/e;
    }
    return %clusters;
}



sub makemalinmatrix () {
    my ($rootsaddress,$malinfile,$infile,$outfile) = @_;
    my @roots = @$rootsaddress;
    my @data = ();
    my $ind;
    $/ = "\n";
    open (IN, $malinfile);
    $_ = <IN>;
    ($positions) = /\s+\d+\s+(\d+)/;


    open (IN, $infile);
    open (OUT, ">$outfile");
    while (<IN>) {
	@c = split;
	for $ind (0..$#c) {
	    $data[$ind] .= int ($c[$i] ne "-");
	}
    }
    for $ind (0..$#roots) {
	$data[$ind] .= "0" x ($positions-length $data[$ind]);
	print OUT "$roots[$ind]\t$data[$ind]\n";
    }
    close OUT;
}
	    
    


sub makeexpobsphase() {
    my ($infile,$pointer,$outfile) = @_;
    my %tetranucphasedist = %$pointer;
    open (IN, $infile);
    local $/ = ">";
    my %obsphase = ();
    my %expphase = ();
    # Temporary
    my $ch = 0;
    while (<IN>) {
	next unless /\n[A-Z]+([ACGT]{2})[a-z]+([ACGT]{2})/;
	$b = $1.$2;
	($phase,$group) = /:ph(\d):.*:(\S+)/;
	$obsphase{$group}[$phase]++;
	for $ph (0,1,2) {
	    $expphase{$group}[$ph] += $tetranucphasedist{$b}[$ph];
	}

    }
    open (OUT, ">$outfile");
    print OUT "GROUP\tObs0\tObs1\tObs2\tExp0\tExp1\tExp2\n";
    for $group (sort keys %obsphase) {
	print OUT "$group\t";
	for $ph (0,1,2) {
	    printf OUT "%.4f\t", $obsphase{$group}[$ph];

	}
	for $ph (0,1,2) {
	    printf OUT "%.4f\t", $expphase{$group}[$ph];
	}
	print OUT "\n";
    }
    close OUT;

	
}

sub maketetranuc () {
    my $infile = $_[0];
    open (IN, $infile);
    local $/ = "\n";
    my %tetranucdist = ();
    #  Temporary
    my $ch = 0;
    while (<IN>) {
	next if />/;
	last if $ch++ > 100;
	$_ = join ("", /[A-Z]+/g);
	while (/(?=([ACGT]{4}))/g) {
	    $tetranucdist{$1}[(2 + length $PREMATCH)%3]++;
	}
    }
    for $tn (keys %tetranucdist) {
	$tot = 0;
	for $ph (0,1,2) {
	    $tot += $tetranucdist{$tn}[$ph];
	}
	for $ph (0,1,2) {
	    $tetranucdist{$tn}[$ph] /= $tot;
	}
    }
    return \%tetranucdist;
}




sub compilesharedintrons() {
    my ($root, $rootsaddress, $suffix, $cladesfile, $outfile) = @_;
    my @roots = @$rootsaddress;
    my $index = "";


    my %divfromroot = ();
    my %sharingcompilation = ();
    open (IN, "$cladesfile") || die;
    local $/ = "\n";
    while (<IN>) {
	/Divergence:\s+$root\s+(\S+)\s+(\d+)/ && 
	    ($divfromroot{$1} = $2);
    }
    close IN;
	
    for $r2 (@roots) {
	next if $r2 eq $root;
	if ( (open (IN, "$root-$r2.sharedintrons"))
	    &&
	     (open (IN, "$r2-$root.sharedintrons"))) {
	    die "### FATAL ERROR\nYou've got both $root-$r2.sharedintrons and $r2-$root.sharedintrons.  Concerned one might be empty?";
	}
	if (open (IN, "$root-$r2.sharedintrons")) {$index = 0;} 
	elsif (open (IN, "$r2-$root.sharedintrons")) {$index = 1; }


	else {
	    die "### FATAL ERROR\nCan't find $root-$r2.sharedintrons or $r2-$root.sharedintrons\n"}
	
	while (<IN>) {
	    @F = split;
	    unless ($F[$index] eq "-") {
		$sharingcompilation{$F[$index]}{$divfromroot{$r2}}{("Y","N")[$F[int !$index] eq "-"] .":$r2:"}++;
	    }
	}
    }
    open (OUT, ">$outfile");
    for $intron (sort keys %sharingcompilation) {
	print OUT "$intron\t";
	for $clade (sort {$a<=>$b} keys %{$sharingcompilation{$intron}}) {
	    print OUT "$clade:(";
	    print OUT join (",", keys %{$sharingcompilation{$intron}{$clade}});
	    print OUT ") ";
	}
	print OUT "\n";
    }
    close OUT;
    
}
	



sub makeintronsflankswithsharing() {
    ($infileflanks,$infileintronsharing,$outfile) = @_;
    #  Note that here we are using parsimony.  If it's in a species
    #  from which our species diverged at a given node, we place it in that node.
    my %pattern = ();
    open (IN, $infileintronsharing);
    local $/ = "\n";
    while (<IN>) {
	
	($int,@sharing) = split;
	for $info (@sharing) {
	    $node = int $info;
	    %found = $info =~ /[(,]([A-Z])(:)/g;
	    if (exists $found{"Y"}) {
		$pattern{$node}{$int} = "Y";
	    }
	    elsif (exists $found{"N"}) {
		$pattern{$node}{$int} = "N";
	    }
	}
    }
    open (IN, $infileflanks);
    open (OUT, ">$outfile");
    local $/ = ">";
    while (<IN>) {
	chomp;
	next unless /\S/;
	($head,$int,$seq) = /((\S+?):\S+).*\n(\S+)/;
	$patt = "";
	$nodeorder = join (" ", 0, sort keys %pattern);
	for $node (sort keys %pattern) {
	    if (exists $pattern{$node}{$int}) {
		$patt .= "$node$pattern{$node}{$int}-";
	    }
	    else { $patt .= "$node?-";}
	}
	$patt =~ s/-$//;
	if ($patt =~ /.*(\d)Y/) {
	    $gainposmin = $1;
	}
	else { $gainposmin = 0;}
	if ($patt =~ /(\d+)N[^Y]*$/) {
	    $not = $1;
	    $nodeorder =~ /(\d+) $not/;
	    $gainposmax = $1;

	}
	else { $gainposmax = 4;}
	print OUT ">$head:$patt:$gainposmin-$gainposmax\n$seq\n";
    }
    close OUT;
	
}

sub makeexonsintrons () {
    my ($infilegff,$infilefasta,$outfile) = @_;
    open (IN, $infilegff);
    local $/ = "\n";
    my %coords = ();
    my %strand = ();
    while (<IN>) {
	next unless /^\S+\s+\S+\s+CDS\s/;
	@F = split;
	#  Get gene name... 
	if (
	    /Parent=([^\s;]+)/ ||
	    /transcript_id "(\S+)"/ ||
	    /transcriptId (\S+)/ ||
	    
	    /mRNA (\S+);/) {
	    $id = $1;
	    $id =~ s/.+\|//;
	    $id =~ s/:/_/g;
	    $coords{$F[2]}{$F[0]}{$id}{"$F[3]..$F[4]"}++;
	    $strand{$id} = $F[6];
	}
	
    }
    if (exists $coords{"CDS"}) {
	%coords = %{$coords{"CDS"}};
    }
    else { %coords = %{$coords{"exon"}}}
    open (IN, "$infilefasta");
    open (OUT, ">$outfile");
    local $/ = ">";
    while (<IN>) {
	($n) = /(\S+).*/;
	chomp;
	$seq = join ("", $POSTMATCH =~ /\S+/g);
	for $g (keys %{$coords{$n}}) {
	    $finalcoords = $coords = join (",", sort {$a<=>$b} keys %{$coords{$n}{$g}});
	    ($l,$r) = ($coords =~ /\d+/g)[0,-1];
	    $geneseq = substr ($seq,$l-1,$r-$l+1);
	    $geneseq =~ tr/A-Z/a-z/;
	    $coords =~ s/\d+/$&-$l/eg;
	    while ($coords =~ /(\d+)\.\.(\d+)/g) {
		$bit = substr ($geneseq,$1,$2-$1+1);
		$bit =~ tr/a-z/A-Z/;
		substr ($geneseq,$1,$2-$1+1,$bit);
	    }
	    if ($strand{$g} eq "-") {
		$geneseq = reverse $geneseq;
		$geneseq =~ tr/ACGTacgt/TGCAtgca/;
	    }
	    print OUT ">$g\t$strand{$g}$n:$finalcoords\n$geneseq\n";
	}
    }
    close OUT;
}


sub checkforfiles() {
    my @roots = @_;
    my $r;
    for $r (@roots) {
	open (IN, "$r.gff") || die "###FATAL ERROR\n$r.gff not found.\nSpeciesName.gff and SpeciesName.fna needed for each species.";
	open (IN, "$r.fna") || die "###FATAL ERROR\n$r.fasta not found.\nSpeciesName.gff and SpeciesName.fna needed for each species.";
    }
}

	
sub maketree() {
    #(1,((2,3),4))
    my ($infile,$outfile) = @_;
    open (IN, $infile);
    local $/ = undef;
    my $tree = <IN>;
    $tree= join (" ", "", $tree =~ /([^,\(\)]+|\(|\))/g, "");
    my ($bit,$b) = ();
    while ($tree =~ /(?<=\s)[()]/g) {
	$bit = $&;
	$b += 2* ($& eq "(") - 1;
	$tree =~ s//sprintf ":$bit:%d", $b + ($bit eq ")")/e;
    }
    my %relationships = ();
    my @clades = ();
    while ($tree =~ /(?=:\((:\d+ )(.*?):\)\1)/g) {
	$clade = $2;
	$clade =~ s/:[()]:\d+//g;
	push (@clades, join (" ", $clade =~ /\S+/g));
    }
    @clades = sort {$a =~ tr/ / / <=> $b =~ tr/ / /} @clades;
    my %alltaxa = ();
    for $i (0..$#clades) {
	@taxa = $clades[$i] =~ /\S+/g;
	for $t1 (@taxa) {
	    $alltaxa{$t1}++;
	    for $t2 (@taxa) {
		last if $t1 eq $t2;
		unless (exists $divergence{$t1}{$t2}) {
		    $divergence{$t1}{$t2} = $divergence{$t2}{$t1} = $i+1;
		}
	    }
	}
    }
    open (OUT, ">$outfile");
    for $i (0..$#clades) {
	printf OUT "Clade %d:\t$clades[$i]\n", $i+1;
    }
    for $t1 (keys %divergence) {
	for $t2 (keys %divergence) {
	    next if $t1 eq $t2;
	    print OUT "Divergence:\t$t1\t$t2\t$divergence{$t1}{$t2}\n";
	}
    }
    close OUT;
    return (\%divergence, \%alltaxa);
}
    

sub makeorthologsmultispecies () {

    ($rootsaddress, $outfile) = @_;
    local $/ = "\n";
    my @r = @$rootsaddress;
    my %orths = ();
    for $r1 (@r) {
	for $r2 (@r) {
	    last if $r1 eq $r2;
	    if (open (IN, "$r1-$r2.orthologs")) {
		@s = ($r1,$r2);
	    }
	    elsif (open (IN, "$r2-$r1.orthologs")) {
		@s = ($r2,$r1);
	    }
	    else { die "Cannot open $r1-$r2.orthologs or $r2-$r1.orthologs"; }
	    while (<IN>) {
		for $i (0,1) {
		    $orths{$s[$i]}{(split)[$i]}{$s[!$i]} = (split)[!$i];
		    $genetospecies{(split)[$i]} = $s[$i];

		};

	    }
	}
    }
    open (OUT, ">$outfile");
    $r1 = $r[0];
    for $g1 (keys %{$orths{$r1}}) {

	$no = 0;
	%slice = %{$orths{$r1}{$g1}};
	next unless (scalar keys %slice) == $#r;
	for $i2 (1..$#r) {
	    $r2 = $r[$i2];
	    for $i3 ($r2+1..$#r) {
		$r2 = $r[$i2];
		$no++ unless $orths{$r2}{$slice{$r2}}{$r3} eq $slice{$r3};
	    }
	}
	print OUT join ("\t", (reverse sort {$genetospecies{$a} cmp $genetospecies{$b}} ($g1, values %slice)), "\n") unless $no;
	    
    }
    close OUT;
}

		


sub makesharedintrons() {
    my ($ei1,$ei2,$pro1,$pro2,$orths,$outfile) = @_;
    my %orths = ();
    my $filevar = "";
    my %pro = ();
    my %ei = ();
    my %get = ();
    
    #return if open (IN,$outfile);
    open (IN, "$orths");
    open (OUT, ">$outfile");
    local $/ = "\n";
    while (<IN>) {
	$orths{(split)[0]} = (split)[1];
	$get{(split)[0]}++;
	$get{(split)[1]}++;
    }

    for $filevar ($pro1,$pro2) {
	local $/ = ">";
	open (IN, $filevar);
	while (<IN>) {
	    print;
	    /(\S+).*\n(\S+)/;
	    $get{$1} && ($pro{$1} = $2);
	}
    }
    for $filevar ($ei1,$ei2) {
	local $/ = ">";
	open (IN, $filevar);
	while (<IN>) {
	    /(\S+).*\n(\S+)/;
	    $get{$1} && ($ei{$1} = $2);
	}
    }

    for $key (keys %orths) {
	@o = ($key,$orths{$key});
	$command =  "echo '>1\n$pro{$o[0]}\n>2\n$pro{$o[1]}\n' | clustalo --wrap 1000000 -i - | perl -ane 'print unless \$. % 2'";
	(@a) = (qx($command));

	@confidentpositions = &alignmentquality(@a);

	%apos = ();
	for $i (0,1) {
	    $temp = $ei{$o[$i]};

	    $temp =~ s/[a-z]+([A-Z])/lc($1)/eg;

	    @codons = $temp =~ /.../g;
	    $pos = 0;
	    $a[$i] =~ s/[A-Z\*]/$codons[$pos++]/eg;
	    $a[$i] =~ s/-/---/g;
	    $inumb = 0;
	    while ($a[$i] =~ /[a-z]/g) {
		$inumb++;
		$apos{sprintf "%d.%d", (length $PREMATCH)/3+1, (length $PREMATCH) % 3}[$i] = $inumb;
	    }
	}

	for $ipos (sort {$a<=>$b} keys %apos) {
	    if ($confidentpositions[$ipos]) {
		for $i (0,1) {
		    print OUT ("-\t","$o[$i].$apos{$ipos}[$i]\t")[$apos{$ipos}[$i] > 0];
		}
		print OUT "\n";
	    }

	}


    }
    close OUT;
}




sub makesharedintronsmanypaml() {
    my ($rootsaddress,$infileorths,$outfile,$pamloutfile) = @_;
    my @roots = @$rootsaddress;
    my %group = ();
    my $totalconfidentpositions = 0;
    print "Here we are\n";
    open (IN, "$infileorths");
    open (OUT, ">$outfile");
    local $/ = "\n";
    while (<IN>) {
	$group++;
	@F = split;
	for $i (0..$#roots) {
	    $group[$group]{$roots[$i]} = $F[$i];
	    $get{$F[$i]}++;
	}
    }
    print "Read files?\n";

    for $r (@roots) {
	open (IN, "$r.pro");
	local $/ = ">";
	while (<IN>) {
	    /(\S+).*\n(\S+)/;
	    $get{$1} && ($pro{$1} = $2);
	}
	open (IN, "$r.exons-introns");
	local $/ = ">";
	while (<IN>) {
	    /(\S+).*\n(\S+)/;
	    $get{$1} && ($ei{$1} = $2);
	}
    }
    print "But can we read protein and exint files?\n";
    for $group (0..$#group) {
        #print "$group\n";
	%apos = ();
	$string = "";
	for $r (@roots) {

	    $string .= ">$r\n$pro{$group[$group]{$r}}\n";
            #print "$string\n";

	}

	$command = "echo '$string' | clustalo --wrap 1000000 -i - | perl -ane 'print unless \$. % 2'";
        print "MADEIT\n";
	(@a) = (qx($command));
	@confidentpositions = &alignmentqualitymany(@a);

	$totalconfidentpositions += sum(@confidentpositions);	
	
	#  Here we are...
	%thispaml = ();
	for $i (0..$#roots) {
	    $r = $roots[$i];
	    $temp = $ei{$group[$group]{$r}};
	    $temp =~ s/[a-z]+([A-Z])/lc($1)/eg;
	    @codons = $temp =~ /.../g;
	    $pos = 0;
	    $a[$i] =~ s/[A-Z\*]/$codons[$pos++]/eg;
	    $a[$i] =~ s/-/---/g;

	    $inumb = 0;

	    while ($a[$i] =~ /[a-z]/g) {
		$inumb++;
		$apos{sprintf "%d.%d", (length $PREMATCH)/3+1, (length $PREMATCH) % 3}[$i] = $inumb;
	    }
	    $a[$i] =~ tr/a-z/A-Z/;
	    @acodons = $a[$i] =~ /.../g;
	    for $pos (0..$#confidentpositions) {
		($thispaml{$r} .= $acodons[$pos]) 
		    if $confidentpositions[$pos];
	    }
	}
	@existsconservedposition = ();
	for $ipos (sort {$a<=>$b} keys %apos) {
	    for $species (0..$#roots) {
		($existsconservedposition[$species]=1) if $apos{$ipos}[$species];
	    }
	}

	next unless (sum(@existsconservedposition) == scalar @roots);

	for $r (keys %thispaml) {
	    $paml{$r} .= $thispaml{$r};
	}
	for $ipos (sort {$a<=>$b} keys %apos) {
	    if ($confidentpositions[$ipos]) {
		#|| 
		#($apos{$ipos}[0] && $apos{$ipos}[1])) {
		for $i (0..$#roots) {
		    print OUT ("-\t","$group[$group]{$roots[$i]}.$apos{$ipos}[$i]\t")[$apos{$ipos}[$i] > 0];
		}
		print OUT "\n";
	    }

	}
	
	
    }
    close OUT;
    open (OUT, ">$pamloutfile");
    printf OUT "      %d  %d\n", (scalar keys %paml), length $paml{$roots[0]};
    for $r (keys %paml) {
	print OUT "$r  $paml{$r}\n";
    }
    close OUT;
    return length $paml{$roots[0]}
}


    
	


	
    
	
    

sub roots () {
    my @roots = ();
    ($suffix,@files) = @_;
    for $f (@files) {
	$f =~ /(\S+)$suffix/;
	push (@roots, $1);
    }

    @roots;
}


sub makedmnd() {
    my $infile = $_[0];
    system "diamond makedb --in $infile --db $infile --quiet";
}

sub makeorthologs() {
    my ($infile1,$infile2,$outfile) = @_;
    system "diamond blastp -e1E-10 --query $infile1 --db $infile2 | sort --k=12nr | perl -ane 'print unless \$a{\$F[0]} || \$a{\$F[1]}; \$a{\$F[0]}++; \$a{\$F[1]}++' | cut -f 1,2  > $outfile";
}
    


sub makeintronsflanks () {
    my $minflanklength = 10; # Parameter
    my $maxflanklength = 20; # Parameter
    my $minintronlength = 30; # Parameter
    my ($infile,$outfile) = @_;
    open (IN, $infile);
    open (OUT, ">$outfile");
    $/ = ">";
    while (<IN>) {
	($n,$seq) = /(\S+).*\n(\S+)/;
	$seq = join (" ", "", $seq =~ /(?:[a-z]*[A-Z]){3}/g);

	$intnumb = 0;
	$seq =~ s/ ([A-Z]*)(?=[a-z]+)/sprintf " $1:%d:%d:", ++$intnumb, int length $1/eg;
	$seq = join ("", $seq =~ /\S+/g);
	while ($seq =~ /([A-Z]{$minflanklength,$maxflanklength})(?=:(\d+):(\d+):([a-z]{$minintronlength,})([A-Z]{$minflanklength,$maxflanklength}))/g) {
	    printf OUT ">$n.$2:%d-%d-%d:ph$3\n$1$4$5\n", length $1,length $4,length $5;
	}
    }
    close OUT;
}

sub makeends () {
    my $maxendlength = 100; #Parameter
    my ($infile, $outfile) = @_;
    local $/ = ">";
    open (IN, $infile);
    open (OUT, ">$outfile");
    while (<IN>) {
	chomp;
	next unless /\S/;
	($n,$e1,$i,$e2) =  /(\S+)\n([A-Z]+)([a-z]+)([A-Z]+)/;
	$i =~ /.{1,$maxendlength}/;
	print OUT ">$n:up\n$e1$&\n";
	$i =~ /.{1,$maxendlength}$/;
	print OUT ">$n:down\n$&$e2\n";
    }
    close OUT;
}
	
sub blastends () {
    my ($infile, $outfile) = @_;
    
    system "makeblastdb -in $infile -dbtype nucl";
    
    system "megablast -m 8 -i $infile -d $infile -e0.00001 > $outfile";
}


sub makepro () {
    my %t = ('TCA','S','TCC','S','TCG','S','TCT','S','TTC','F','TTT','F','TTA','L','TTG','L','TAC','Y','TAT','Y','TAA','*','TAG','*','TGC','C','TGT','C','TGA','*','TGG','W','CTA','L','CTC','L','CTG','L','CTT','L','CCA','P','CCC','P','CCG','P','CCT','P','CAC','H','CAT','H','CAA','Q','CAG','Q','CGA','R','CGC','R','CGG','R','CGT','R','ATA','I','ATC','I','ATT','I','ATG','M','ACA','T','ACC','T','ACG','T','ACT','T','AAC','N','AAT','N','AAA','K','AAG','K','AGC','S','AGT','S','AGA','R','AGG','R','GTA','V','GTC','V','GTG','V','GTT','V','GCA','A','GCC','A','GCG','A','GCT','A','GAC','D','GAT','D','GAA','E','GAG','E','GGA','G','GGC','G','GGG','G','GGT','G',);
    

    my ($infile, $outfile) = @_;
    open (IN, $infile);
    open (OUT, ">$outfile");
    local $/ = ">";
    while (<IN>) {
	chomp;
	next unless /\S/;
	($n,$seq) = /(\S+).*\n(\S+)/;

	$pos = 0;
	$intpos = "";
	@e = $seq =~ /[A-Z]+/g;

	$seq = join ("", $seq =~ /[A-Z]/g);
	$pro = "";
	while ($seq =~ /[A-Z]{3}/g) {
	    if (exists $t{$&}) {$pro .= $t{$&}}
	    else {$pro .= "X"}
	}
	print OUT ">$n\n$pro\n";
    }
    close OUT;
}


sub makeparalogs () {
    my ($infile, $outfile) = @_;
    system "diamond makedb --in $infile --db $infile";
    system "diamond blastp -e1E-20 --query $infile --db $infile | cut -f 1,2 > $outfile";
}


sub filterblastends () {
    my ($infile, $paralogfile, $outfile) = @_;
    my %p = ();
    open (IN, $paralogfile);
    local $/ = "\n";
    while (<IN>) {
	$p{(split)[0]}{(split)[1]}++;
	$p{(split)[1]}{(split)[0]}++;
    }
    open (IN, $infile);
    local $/ = "\n";
    open (OUT, ">$outfile");
    while (<IN>) {
	@g = /\b(\S+?)\.\d+:/g;
	@i = /\b(\S+?):/g;
	($i[0] eq $i[1]) || 
	    $p{$g[0]}{$g[1]} || print OUT $_;
    }
    close OUT;
}


sub makematches () {
    print "Hello\n";
    my ($infile, $outfile) = @_;
    open (IN, $infile);
    local $/ = "\n";
    open (OUT, ">$outfile");
    while (<IN>) {
	@F = split;
	($qn,$qu,$qi,$qd,$qb) = $F[0] =~ /(\S+?:(\d+)-(\d+)-(\d+):\S+):(up|down)/;
	$qi = &min ($qi,$maxlength);
	($sn,$su,$si,$sd,$sb) = $F[1] =~ /(\S+?:(\d+)-(\d+)-(\d+):\S+):(up|down)/;
	$si = &min ($si,$maxlength);
	$sameend = int ($qb eq $sb);
	$samestrand = int ($F[9]>$F[8]);
	next if (($sameend + $samestrand)%2);
	if ($qb eq "up") {
	    next unless $F[6] < $qu+10;
	    next unless $F[6] > $qu-5;
	}
	else {
	    next unless $F[7] > $qi-10;
	    next unless $F[7] < $qi+5;
	}
	
	if ($sb eq "up") {
	    $near = &min(@F[8,9]);
	    next unless $near < $su+10;
	    next unless $near > $su-5;
	}
	else {
	    $near = &max(@F[8,9]);
	    next unless $near > $si-10;
	    next unless $near < $si+5;
	}
	
	$strand = ("-","+")[$sb eq $qb];
	
	if ($qb eq "down") {
	    
	    if (exists $matches{$qn}{$sn}{$strand}) {
		print OUT "$qn\t$sn\t$strand\n";
	    }
	}
	else {
	    $matches{$qn}{$sn}{$strand}++;
	}
    }
    close OUT;
}


    
sub makeclusters () {		  
    my ($infile, $prefix) = @_;
    my %m = ();
    open (IN, $infile);
    local $/ = "\n";
    while (<IN>) {
	 /(\S+?)\.\d+:\S+\s+(\S+?)\.\d+:\S+\s+/;
	 $m{(split)[0]}{(split)[1]} = (split)[2];
	 $m{(split)[1]}{(split)[0]} = (split)[2];
    }
    while (scalar keys %m) {
	$k1 = (keys %m)[0];
	do {
	    $new = 0;
	    for $k2 (keys %{$m{$k1}}) {
		for $k3 (keys %{$m{$k2}}) {
		    next if exists $m{$k1}{$k3};
		    $m{$k1}{$k3} = ("-","+")[int ($m{$k1}{$k2} eq $m{$k2}{$k3})];
		    $new++;
		}
	    }
	} while ($new);
	$group++;
	$vals = join ("", values %{$m{$k1}});
	$neg = ("+","-")[int ($vals =~ tr/-/-/)  > ($vals =~ tr/+/+/)];
	$group{$k1} = $neg . $group;
	for $k2 (keys %{$m{$k1}}) {
	    next if $k1 eq $k2;
	    $str = ("-","+")[int ($m{$k1}{$k2} eq $neg)];
	    $group{$k2} = $str.$group;
	    delete $m{$k2};
	}
	delete $m{$k1};
    }
    return %group;
}

sub printclusters () {
    my ($infile, $outfile1, $outfile2, %group) = @_;

    
    open (IN, $infile);
    open (OUT, ">$outfile1");
    local $/ = ">";
    my %printing = ();
    while (<IN>) {
	chomp;
	($n) = /(\S+)/;
	if (exists $group{$n}) {
	    
	    s/(?<=\S)(?=\s)/:Group$group{$n}/;
	    $printing{$group{$n}} .= ">$_";
	}
	else {
	    s/(?<=\S)(?=\s)/:NoGroup/;
	}
	print OUT ">$_";
    }
    close OUT;
    open (OUT, ">$outfile2");
    print join (" ", %printing);
    for $group (sort {abs($a)<=>abs($b)} keys %printing) {
	print OUT "## Group $group\n$printing{$group}\n";
    }
    close OUT;
}


sub alignmentquality () {
    my $window = 10;  # Parameter
    my $threshold = 0.4; # Parameter
    my @confidentpositions = ();


    my @a1 = $_[0] =~ /\S/g;
    my @a2 = $_[1] =~ /\S/g;

    $min = int $window*$threshold;
    for $i (0..$#a1) {
	my ($matchesleft,$matchesright) = (0,0);
	for $j (max(0,$i-$window)..$i-1) {
	    $matchesleft += $a1[$j] eq $a2[$j];
	}
	for $j ($i+1..$i+$window) {
	    $matchesright += $a1[$j] eq $a2[$j];
	}
	$confidentpositions[$i] = int (
	    ($matchesleft >= $threshold*$window) &&
	    ($matchesright >= $threshold*$window));
    }

    return (@confidentpositions[0..$#confidentpositions]);
}    
    

    
    


sub alignmentqualitymany() {
    my @a = @_;

    @qcompile = ();
    for $i1 (0..$#a) {
	for $i2 ($i1+1..$#a) {

	    @q = &alignmentquality(@a[$i1,$i2]);
	    @qcompile = &addarrays(\@q,\@qcompile);
	}
    }
    for $index (0..$#qcompile) {
	$qcompile[$index] = ($qcompile[$index] == $#a*($#a+1)/2);
    }


    @qcompile;
}

sub addarrays() {
    my ($a1a, $a2a) = @_;
    my @a1 = @$a1a;
    my @a2 = @$a2a;
    my @a3 = ();

    for $i (0..$#a1) {
	$a3[$i] = $a1[$i]+$a2[$i];
    }
    @a3;
}
	





sub min () {
    (sort {$a<=>$b} @_)[0];
}


sub max () {
    (sort {$a<=>$b} @_)[-1];
}

	  

sub sum () {
    my $sum = 0;
    my $bit;
    for $bit (@_) {
	$sum += $bit;
    }
    return $sum;
}

	

