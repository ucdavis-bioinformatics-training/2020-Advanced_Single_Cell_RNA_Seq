#!/usr/bin/perl

$inb=0;
$rb=0;
open ($fp, "<$ARGV[0]");
while ($line=<$fp>) {
    chomp $line;

    if ($line =~ /^```\{r/) {
        $rb = 1;
    }

    if ($line =~ /^```$/) {
        if ($inb==0 && $rb==0) {
            print "<div class=\"output\">";
            $inb=1;
            next;
        } elsif ($inb==0 && $rb==1) {
            $rb=0;
        } else {
            print "</div>\n";
            $inb=0;
            next;
        }
    }

    if ($inb==1) {
        $line=~s/^## //;
    }
    print "$line\n";
}
close($fp);
