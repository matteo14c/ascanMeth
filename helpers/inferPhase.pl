open(OUT,">phase.err.log");
$vcf=shift;
open(IN,$vcf);
while(<IN>)
{
	if ($_=~/^#/)
	{
		print;
		next;
	}

	($chr,$start,$p,$ref,$alt,$score,$filter,$f,$annots_F,$Cd,$P1,$P2)=(split());
	next unless $score>=20;
	($C,@Can)=(split(/\:/,$Cd));
	($P1,@P1an)=(split(/\:/,$P1));
	($P2,@P2an)=(split(/\:/,$P2));
	next unless $C=~/\|/;
	($c1,$c2)=   $C =~ /\|/ ? (split(/\|/,$C)) : (split(/\//,$C)) ;
	($p11,$p12)= $P1 =~ /\|/ ? (split(/\|/,$P1)) : (split(/\//,$P1)) ;
	($p21,$p22)= $P2 =~ /\|/ ? (split(/\|/,$P2)) : (split(/\//,$P2)) ;

	$Can=join(":",@Can);
	$P1an=join(":",@P1an);
	$P2an=join(":",@P2an);	
	
	$outH1=0;
	$outH2=0;
	$sep="|";

	$c1P1= $c1 == $p11 || $c1 == $p12 ? 1 : 0;
	$c2P1= $c2 == $p11 || $c2 == $p12 ? 1 : 0;
	$c1P2= $c1 == $p21 || $c1 == $p22 ? 1 : 0;
	$c2P2= $c2 == $p21 || $c2 == $p22 ? 1 : 0;


	if ($c1P1==1 && $c1P2==0)
	{
		$outH1=$c1;
		$outH2=$c2;
		$phased=1
	}elsif ($c1P1==0 && $c1P2==1){
		$outH2=$c1;
		$outH1=$c2;
		$phased=1;
	}elsif ($c2P1==1 && $c2P2==0){
		$outH1=$c2;
                $outH2=$c1;
                $phased=1
	}elsif($c2P1==0 && $c2P2==1){
		$outH2=$c2;
                $outH1=$c1;
                $phased=1
	}else{
		$outH1=$c1;
		$outH2=$c2;
		$sep="/";
	}
	$out="$outH1$sep$outH2";
	print OUT "$chr\t$start\t$end\t$C\t$P1\t$P2\t$out:inv\n" if $out ne $C;
	print OUT "$chr\t$start\t$end\t$C\t$P1\t$P2\t$out:cor\n" if $out eq $C;
	print "$chr\t$start\t$p\t$ref\t$alt\t$score\t$filter\t$f\t$annots_F\t$out:$Can\t$P1:$P1an\t$P2:$P2an\n";
}
