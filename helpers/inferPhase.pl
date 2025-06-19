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
	#next unless $C=~/\|/;
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
		if ($c1 eq $p11)
		{
			$pH1=$p11;
			$pH2=$p12;
		}else{
			$pH1=$p12;
			$pH2=$p11;
		}

	}elsif ($c1P1==0 && $c1P2==1){
		$cH2=$c1;
		$cH1=$c2;
		if ($c1 eq $p21)
                {
                        $mH1=$p22;
                        $mH2=$p21;
                }else{
			$mH1=$p21;
                        $mH2=$p22;
                }
	
	
	}elsif ($c2P1==1 && $c2P2==0){
		$cH1=$c2;
                $cH2=$c1;

		if ($c2 eq $p11)
                {
                        $pH1=$p11;
                        $pH2=$p12;
                }else{
                        $pH1=$p12;
                        $pH2=$p11;
                }



	}elsif($c2P1==0 && $c2P2==1){
		$cH2=$c2;
                $cH1=$c1;

		if ($c2 eq $p21)
                {
                        $mH1=$p22;
                        $mH2=$p21;
                }else{
                        $mH1=$p21;
                        $mH2=$p22;
                }


	}else{
		$cH1=$c1;
		$cH2=$c2;
		$sep="/";
	}
	$C_m="$cH1$sep$cH2";
	$P_m="$pH1$sep$pH2";
	$M_m="$mH1$sep$mH2";
	
	$tag_er="";
	$tag_er.="_inv" if $C_m ne $C && $sep eq "|";
	$tag_er.="_noP" if $c1P1 ==0 && $c2P1==0; # nessun allele paterno
        $tag_er.="_noM" if $c1P2 ==0 && $c2P2==0; # nessun allele paterno

	#quick note:
	# right now I miss out the cases where C is homozygous and the parents are hets
	# this could be phased in the parents
	# quick fix here
	if ($C_m eq "1/1" && ($P_m eq "0/1" || $P_m eq "1/0") && ($M_m eq "0/1" || $M_m eq "1/0") )
	{
		$sep="|";
		$P_m="1|0";
		$M_m="0|1";

	}elsif( $C_m eq "0/0" && ($P_m eq "0/1" || $P_m eq "1/0") && ($M_m eq "0/1" || $M_m eq "1/0") ){
		$sep="|";
                $P_m="0|1";
                $M_m="1|0";

	}

	if ($tag_er=~/_noM/)
	{
		#what to do?
	}elsif($tag_er=~/_noP/){
		#what to do?
	}
	print OUT "$chr\t$start\t$end\t$C\t$P1\t$P2\t$C_m:$tag_er\n" if $tag_er ne "";
	#$P1an, $P2an,$Can
	print "$chr\t$start\t$p\t$ref\t$alt\t$score\t$filter\t$f\t$annots_F\t$C_m\t$P_m\t$M_m\n" if $sep eq "|";
}
