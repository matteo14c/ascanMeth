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
	
	# for child phasing
	$cH1=0;
	$cH2=0;
	$sep="|";

	# for dad's phasing
	$pH1=$p11;
	$pH2=$p12;
	$sep="|";

	# for mom's phasing
	$mH1=$p21;
	$mH2=$p22;

	$c1P1= $c1 == $p11 || $c1 == $p12 ? 1 : 0;
	$c2P1= $c2 == $p11 || $c2 == $p12 ? 1 : 0;
	$c1P2= $c1 == $p21 || $c1 == $p22 ? 1 : 0;
	$c2P2= $c2 == $p21 || $c2 == $p22 ? 1 : 0;


	if ($c1P1==1 && $c1P2==0)
	{
		$cH1=$c1;
		$cH2=$c2;
		$pH1=$p11;
		$pH2=$p12;

		$phased=1
	}elsif ($c1P1==0 && $c1P2==1){
		$cH2=$c1;
		$cH1=$c2;
		$mH1=$p22;
		$mH2=$p21;

		$phased=1;
	}elsif ($c2P1==1 && $c2P2==0){
		$cH1=$c2;
                $cH2=$c1;
                $phased=1;
		$pH1=$p12;
                $pH2=$p11;


	}elsif($c2P1==0 && $c2P2==1){
		$cH2=$c2;
                $cH1=$c1;
                $phased=1;
		$mH1=$p21;
                $mH2=$p22;

	}else{
		$cH1=$c1;
		$cH2=$c2;
		$sep="/";
	}
	$C="$cH1$sep$cH2";
	$P="$pH1$sep$pH2";
	$M="$mH1$seo$mH2";
	
	print OUT "$chr\t$start\t$end\t$C\t$P1\t$P2\t$out:inv\n" if $out ne $C;
	print OUT "$chr\t$start\t$end\t$C\t$P1\t$P2\t$out:cor\n" if $out eq $C;
	#$P1an, $P2an,$Can
	print "$chr\t$start\t$p\t$ref\t$alt\t$score\t$filter\t$f\t$annots_F\t$C\t$P1\t$P2\n";
}
