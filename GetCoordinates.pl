#  Command: perl GetCoords.pl *Pass

@files = @ARGV;

#  Figure out which species each input file is from
for $f (@ARGV) {
    ($r) = $f =~ /(\S+)\.Group\d+\./;
    $files{$r}{$f}++;
}

#  For each species...
for $species (keys %files) {
    #  Open up the corresponding exons-introns file
    open (IN, "$species.exons-introns");
    local $/ = ">";
    %coords = ();
    #  Get the intron coordinates for each intron
    while (<IN>) {
	($g,$s,$contig,$coords) = /(\S+)\s+([\+\-])(\S+):\d+\.\.([\d,\.]+)\.\.\d+/;
$coords =~ s/(\d+),(\d+)/sprintf " %d-%d ",$1+1,$2-1/eg;
@ints = $coords =~ /\d+\-\d+/g;
	if ($s eq "-") { @ints = reverse @ints }
	for $int (0..$#ints) {
	    $coords{sprintf "$g.%d", $int+1} = "$contig:$s$ints[$int]";
	}
    }
    #  For each input file, go through and make a corresponding file that has the
    #  coordinates added.
    for $f (keys %{$files{$species}}) {
	local $/ = ">";
	open (IN, "$f");
	open (OUT, ">$f.withcoords");
	<IN>;
	while (<IN>) {
	    chomp;
	    s/(\S+?):\d+\-\d+\-\d+:\S+/$& $coords{$1}/;
	    print OUT ">$_";
	}
    }
}
