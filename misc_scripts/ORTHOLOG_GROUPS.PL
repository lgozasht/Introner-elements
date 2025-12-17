#!/usr/bin/perl
while (<>) {
    $orth_pairs{join (" ", sort /\S+/g)}++;
}

foreach $pair (keys %orth_pairs) {
    @orths = $pair =~ /\S+/g;
    $orths{$orths[0]} .= "$orths[1] ";
    $orths{$orths[1]} .= "$orths[0] ";
}


foreach $pair (keys %orth_pairs)  {
    &groups($pair);
}

foreach $key (keys %obsolete) {
    delete $groups{$key};
}
print join ("\n", keys %groups);
print "\n";




sub groups () {
    my $current_group = $_[0];
    $groups{$current_group}++ && return;
    foreach $new_memb (&new_membs($current_group)) {

        @members = sort ($new_memb, $current_group =~ /\S+/g);
        for my $i (0..$#members) {
            $obsolete{join (" ", @members[0..$i-1,$i+1..$#members])}++;
        }
        &groups(join (" ", @members));
    }
}
    
sub new_membs () {
    my $orths = my $current_group = $_[0];
    $numb_members = (() = $current_group =~ /\S+/g);

    $orths =~ s/\S+/$orths{$&}/g;
    my %count = ();
    my %success = ();
    foreach $b ($orths =~ /\S+/g) {
	(++$count{$b} == $numb_members) && $success{$b}++;
    }
    keys %success;
}
