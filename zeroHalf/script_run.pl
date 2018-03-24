use strict;
use warnings;
use 5.010;

my @names =( "j12051_2.sm", "j1206_8.sm", "j3014_5.mm", "j3017_4.mm", "j3022_4.mm", "j3029_4.mm", "j309_10.mm", "j309_5.mm", "j6013_3.sm", "j6029_7.sm", "j6045_6.sm", "j6047_3.sm", "j609_1.sm", "j609_3.sm", "j609_8.sm", "j9013_2.sm", "j9021_1.sm", "j9025_2.sm", "j9030_9.sm", "j9041_4.sm", "j909_8.sm");

my $n_1 = 1000; 
my $n_2 = 0;
my $n_3 = 3000;
my $n_4 = 0;
my $n_5 = 0;
my $conc = "";
my $cont=0;
foreach my $n (@names){
       for($n_3 = 100000;$n_3<1000000;$n_3+=100000){
	for($cont = 0; $cont<=15 ;$cont++){
	  for($n_2 = 4; $n_2<=16; $n_2*=2){
	    for($n_4 = 10 ; $n_4 <= 100; $n_4 += 90){	
		$conc = $n."_s";
		$conc .= $cont;
		$conc .=".txt";
		if($n_2 == 4){
		   $n_5 = 3;
		   system("./m2 $conc $n_1 $n_2 $n_3 $n_4 1 $n_5 >>Result/Resultado_$conc"."_$n_1"."_$n_2"."_$n_3"."_$n_4"."_$n_5.txt");	
		   $n_5 = 2;
		   system("./m2 $conc $n_1 $n_2 $n_3 $n_4 1 $n_5 >>Result/Resultado_$conc"."_$n_1"."_$n_2"."_$n_3"."_$n_4"."_$n_5.txt");	
		} 
		if($n_2 == 8){
		   $n_5 = 4;
		   system("./m2 $conc $n_1 $n_2 $n_3 $n_4 1 $n_5 >>Result/Resultado_$conc"."_$n_1"."_$n_2"."_$n_3"."_$n_4"."_$n_5.txt");		
		   $n_5 = 6;
		   system("./m2 $conc $n_1 $n_2 $n_3 $n_4 1 $n_5 >>Result/Resultado_$conc"."_$n_1"."_$n_2"."_$n_3"."_$n_4"."_$n_5.txt");	
		}
		 
		if($n_2 == 16){
		   $n_5 = 10;
		   system("./m2 $conc $n_1 $n_2 $n_3 $n_4 1 $n_5 >>Result/Resultado_$conc"."_$n_1"."_$n_2"."_$n_3"."_$n_4"."_$n_5.txt");		
		   $n_5 = 12;
		   system("./m2 $conc $n_1 $n_2 $n_3 $n_4 1 $n_5 >>Result/Resultado_$conc"."_$n_1"."_$n_2"."_$n_3"."_$n_4"."_$n_5.txt");	
		}
	    }
	  }
	}
	}
}

exit;

