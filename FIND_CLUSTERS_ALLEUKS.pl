#  usage
#  perl FIND_CLUSTERS_ALLEUKS.PL *.gff
#  For each species, there must exist two files
#  SOMETHING.gff -- annotation file, flexibility about specific format
#  SOMETHING.fna -- matching genome fasta file.


use English;
$maxlength = 100;
$minimumfamilies = 4;

#  Section1 : Get relationships from species tree (first argument)

@GFFS = @ARGV;
@roots = &roots(@GFFS);


&checkforfiles(@roots);

foreach $root (@roots) {
    &makeexonsintrons ("$root.gff","$root.fna","$root.exons-introns"); # exons-introns file in $root.exons-introns
    &makeintronsflanks("$root.exons-introns","$root.intronsflanks"); # introns+flanking sequence file in $root.intronsflanks

    &makepro("$root.exons-introns","$root.pro"); # protein file in $root.pro
    &makedmnd("$root.pro"); # make diamond database file for $root.pro
    &makeends("$root.intronsflanks","$root.ends");
    &blastends("$root.ends","$root.ends.blastout");
    &makeparalogs("$root.pro","$root.paralogs");
    &filterblastends("$root.ends.blastout","$root.paralogs", "$root.ends.blastout.filtered");
    &makematches ("$root.ends.blastout.filtered","$root.intronmatches");
    %clusters = &makeclusters("$root.intronmatches");
    %paralogclusters = &makeclusters("$root.paralogs");
    &printclusters ("$root.intronsflanks", $root, \%clusters, \%paralogclusters);

}


	
sub roots () {
    my @roots = ();
    my $f;
    for $f (@_) {
	$f =~ /\S+(?=\.gff)/;
	push (@roots, $&);
    }
    @roots;
}


sub makeexonsintrons () {
    my ($infilegff,$infilefasta,$outfile) = @_;
    open (IN, $infilegff);
    my %coords = ();
    my %strand = ();
    local $/ = "\n";
    while (<IN>) {
	next unless /^\S+\s+\S+\s+(exon|CDS)\s/;
	@F = split;
	#  Get gene name... 
	if (/CDS (\S+)/ ||
	    /Parent=([^\s;]+)/ ||
	    /transcript_id "(\S+)"/ ||
	    /transcriptId (\S+)/ ||
	    /ID=([^\s;]+)/ ||
	    /mRNA (\S+);/) {

	    $coords{$F[2]}{$F[0]}{$1}{"$F[3]..$F[4]"}++;
	    $strand{$1} = $F[6];
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
}


sub checkforfiles() {
    my @roots = @_;
    my $r;
    for $r (@roots) {
	open (IN, "$r.gff") || die "###FATAL ERROR\n$r.gff not found.\nSpeciesName.gff and SpeciesName.fna needed for each species.";
	open (IN, "$r.fna") || die "###FATAL ERROR\n$r.fasta not found.\nSpeciesName.gff and SpeciesName.fna needed for each species.";
    }
}

	

sub addarrays() {
    my ($a1a, $a2a) = @_;
    my @a1 = @$a1a;
    my @a2 = @$a2a;
    my @a3 = ();
    for $i (0..$#a1) {
	$a3[$i] = $a1[$i]+$a2[$i];
    }
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
    local $/ = ">";
    while (<IN>) {
	($n,$seq) = /(\S+).*\n(\S+)/;
	$seq = join (" ", "", $seq =~ /(?:[a-z]*[A-Z]){3}/g);

	$intnumb = 0;
	$seq =~ s/ ([A-Z]*)(?=[a-z]+)/sprintf " $1:%d:%d:", ++$intnumb, int length $1/eg;
	$seq = join ("", $seq =~ /\S+/g);
	while ($seq =~ /([A-Z]{$minflanklength,$maxflanklength}):(\d+):(\d+):([a-z]{$minintronlength,})([A-Z]{$minflanklength,$maxflanklength})/g) {
	    printf OUT ">$n.$2:%d-%d-%d:ph$3\n$1$4$5\n", length $1,length $4,length $5;
	}
    }
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
    close IN;
    open (IN, $infile);
    local $/ = "\n";
    open (OUT, ">$outfile");
    while (<IN>) {
	@g = /\b(\S+?)\.\d+:/g;
	$p{$g[0]}{$g[1]} || print OUT $_;
    }
    close OUT;
    close IN;
}


sub makematches () {
    my ($infile, $outfile) = @_;
    my %matches = ();
    open (IN, $infile);
    open (OUT, ">$outfile");
    local $/ = "\n";
    while (<IN>) {
	@F = split;
	($qn,$qu,$qi,$qd,$qb) = $F[0] =~ /(\S+:(\d+)-(\d+)-(\d+):\S+):(up|down)/;

	$qi = &min ($qi,$maxlength);
	($sn,$su,$si,$sd,$sb) = $F[1] =~ /(\S+:(\d+)-(\d+)-(\d+):\S+):(up|down)/;
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
    my ($infile, $root, $grouppoint, $paralogpoint) = @_;
    my %group = %$grouppoint;
    my %paralogclusters = %$paralogpoint;
    my @families = ();
    system "rm $root.*Group*";
    open (IN, $infile);
    open (OUT, ">$outfile");
    local $/ = ">";
    while (<IN>) {
	chomp;
	($n) = /(\S+)/;
	($gene) = $n =~ /(\S+?)\.\d+:/;
	if (exists $group{$n}) {
	    ($gn) = $group{$n} =~ /(\d+)/g;
	    open (OUT, ">>$root.Group$gn");
	    $paralogclusters{$gene} =~ s/[\+\-]//;
	    $families[$gn]{$paralogclusters{$gene}}++;
	    s/(?<=\S)(?=\s)/:Group$group{$n}:ParalogCluster$paralogclusters{$gene}/;
	}
	else {
	    open (OUT, ">>$root.NoGroup");
	    s/(?<=\S)(?=\s)/:NoGroup/}
	print OUT ">$_";
	close OUT;
    }
    for $group (1..$#families) {
	$numfamilies = scalar keys %{$families[$group]};
	$pass = ("Fail","Pass")[int $numfamilies >= $minimumfamilies];
	system "mv $root.Group$group $root.Group$group.$numfamilies\ParalogFamilies.$pass";
    }
}


sub min () {
    (sort {$a<=>$b} @_)[0];
}


sub max () {
    (sort {$a<=>$b} @_)[-1];
}

	    


