#!/usr/bin/perl
#use warnings;
use strict;

use Bio::SearchIO;
use Getopt::Long;
use Bio::SeqIO;
use Bio::Seq;
use Bio::DB::Fasta;
use Parallel::ForkManager;
 
my $maxProcs = 5;
my $pl = new Parallel::ForkManager($maxProcs);


my (  $input_file 	,
      $output_file 	,
      $frag_min 	,
      $frag_max 	, 
      $rep 		, 
      $repList		,
      $gaps 		, 
      $num_threads 	, 
      $usage 		,
      $overlap 		, 
      $degenerate 	, 
      $cont		,
      $typeModeRun	,
      $mode);
      
      ## pattern values
      $output_file 	= "results";
      $frag_min 	= 2;
      $frag_max 	= 5; 
      $rep 		= 5; 
      $gaps 		= 0;  
      $overlap 		= 0;
      $degenerate 	= 0;
      $mode 		= "n";

GetOptions ( 	"q=s" => \$input_file	, #arquivo fasta
		"o=s" => \$output_file	, #arquivo com as repetitivas--query
                "i=f" => \$frag_min     , #tamanho da subsequencia
                "y=f" => \$frag_max     , #tamanho da subsequencia
		"r=f" => \$rep		, #qtidade minima de repetitivas
		"rl=s"=> \$repList	, #qtidade minima de repetitivas para cada tamanho de subsequencia
		"g=f" => \$gaps		, #qtidade de gaps maximos permitidos entre cada subsequencia
                "t=i" => \$num_threads	, #numero de threads
                "v=f" => \$overlap	, #taxa de overlap permitido
                "d=f" => \$degenerate	, #taxa de degeneraçao
                "m=s" => \$mode		, #rum mode, n to nucleotideo  ou p to protein
		"h"   => \$usage	) ;

if (!$num_threads) { $num_threads = 1 }

if ( defined( $usage ) ) {
    system "clear";
    print " help\n";
    print_usage();
    exit(0);
}


my @repLista = split(/-/,$repList);

##Testa se algum arquivo foi enviado
##se não, apresenta o help e finaliza o processo
if ($input_file eq "" )
{
    system "clear";
    print "NO help\n";
    
    print_usage();
    exit(1);
}

system "clear";
print "RepeatFinder\t Start Process...\n";

##testa se o modo de execução é para proteinas
if($mode eq 'p'){
    $typeModeRun = "repeatFinderProtein";
}else{##caso default DNA
    $typeModeRun = "repeatFinderDNA";
}

for($cont = $frag_min; $cont <= $frag_max; $cont++)
  {
    $pl->start and next;
    my $tmp = $output_file."_".$cont;
    
    ##verifica se a lista de valores tamanho de subsequencia foi declarada e se qtidade é igual a quantidade de subsequencias
    if(defined($repList) && ($#repLista ==($frag_max-$frag_min) )){
	print "./$typeModeRun -q $input_file -o $tmp -i $cont -t $repLista[$cont - $frag_min] -g $gaps -d $degenerate\n";
	system "./$typeModeRun -q $input_file -o $tmp -i $cont -t $repLista[$cont - $frag_min] -g $gaps -d $degenerate";     
    }
    ##verifica se a lista de valores tamanho de subsequencia foi declarada e se qtidade é menor a quantidade de subsequencias
    ##neste caso utiliza somente o primeiro
    elsif(defined($repList) && ($#repLista <($frag_max-$frag_min) )){
	print "./$typeModeRun -q $input_file -o $tmp -i $cont -t $repLista[0] -g $gaps -d $degenerate\n";
	system "./$typeModeRun -q $input_file -o $tmp -i $cont -t $repLista[0] -g $gaps -d $degenerate";     
    }
    ##caso default
    else{
	print "./$typeModeRun -q $input_file -o $tmp -i $cont -t $rep -g $gaps -d $degenerate\n";
	system "./$typeModeRun -q $input_file -o $tmp -i $cont -t $rep -g $gaps -d $degenerate";     
    }
    
    $pl->finish; ## end point of the parallel process
  }
    ##Espera finalizar todos as threads para dar continuidade    
    $pl->wait_all_children;
    
    
##Junto todos os arquivos com resultados em um arquivo unico
system "cat $output_file* > $output_file";


my $output_sts = "Out_$output_file";
open(OUTFILE, ">$output_sts");

my $db = Bio::DB::Fasta->new("$input_file");


exit(0);

#pegar somente a id das sequencias, sem repetição
my @cat_idSequence = `cat $output_file | awk '{printf  \"\%s\\n\", \$1}' | sort -u`;

my $qt = 0;
#exit(0);
print "RepeatFinder\t Processing overlaps..\n";
for(my $cont1 = 1; $cont1 <=$#cat_idSequence; $cont1++){
  
  print "\rRepeatFinder\t Processing : $cont1  of $#cat_idSequence Sequence(s)";
  chomp $cat_idSequence[$cont1];
 
  #faz uma chamada ao S.O. e armazena o resultado em uma hash 
#  print "cat $output_file | awk ' /$cat_idSequence[$cont1]/{printf  \"\%s\\n\", \$0}' | sort -nk5\n";
##  my @cat_result = `cat $output_file | awk ' /$cat_idSequence[$cont1]/{printf  \"\%s\\n\", \$0}' | sort -nk5 `;
  my @cat_result = `cat $output_file | awk '{printf  \"\%s\\n\", \$0}' | grep '$cat_idSequence[$cont1]' |sort -nk5 `;
  
#  cat /usr/lib/cgi-bin/repeatfinder/Out_20131211060941.sts | awk '{print "%s \n", $0}' | grep 'sp|Q8N693|ESX1_HUMAN' | sort -nk5
  
  my (@atual,@proximo);
  my @result;
    
  my $countSequence = $#cat_result;
  #percorre todos os elementos da hash, os quais são cada linha do arquivo output_file
  for (my $contAtual = 0; $contAtual <= $#cat_result; $contAtual++) 
  {
    chomp $cat_result[$contAtual];
     
    #divide a hash em elementos separados e os armazena em um vetor, onde cada um pode ser acessado separadamente
    @atual = split(/[\t]/, $cat_result[$contAtual]); 
    if($atual[5] == ""){next;}
    
    ## se for o ultimo elemento e ele for diferente de nulo,
    ## então insere ele no grupo de respostas.
    if($contAtual == $#cat_result && $atual[5] != ""){
    
	## guarda a sequencia atual como resultado
	$result[$qt] = $cat_result[$contAtual];
	## incrementa o contador de qtidade de microssatelites válidos
	$qt++;
	
	
    }
    
       
    for(my $contProximo = $contAtual + 1 ; $contProximo <= $#cat_result; $contProximo++){
       chomp $cat_result[$contProximo];
       #divide a hash em elementos separados e os armazena em um vetor, onde cada um pode ser acessado separadamente
       @proximo = split(/[\t]/, $cat_result[$contProximo]); 
       
       if($proximo[4] == ""){next; };
      
      if($proximo[4] < $atual[5]){
	#existe um overlap
	
	my $aux_tamanho = $atual[5] - $proximo[5];
	my $overlap_tmp;
	##Verifica se o proximo tem end position maior ou menor que o atual
	if( $aux_tamanho > 0){
	  ## calcula o valor do overlap/ no caso que o proximo tem posição final menor que a posição final do atual
	  $overlap_tmp = int((($atual[5] - $proximo[4] - $aux_tamanho)/($atual[5] - $atual[4]))*100);
	}else{
	  ## calcula o valor do overlap/ no caso que o proximo tem posição final maior que a posição final do atual
	  $overlap_tmp = int((($atual[5] - $proximo[4])/($atual[5] - $atual[4]))*100);
	}
	

	## verifica se o microssatelite atual é menor que o proximo
	## se sim, calcula o overlap do atual tomando como referencia o tamanho do proeximo(isto pq ele é o maior)
	if(($atual[5] - $atual[4]) < ($proximo[5] - $proximo[4])){
		## calcula o valor do overlap/ no caso que o proximo tem posição final maior que a posição final do atual
		$overlap_tmp = int((($atual[5] - $proximo[4])/($proximo[5] - $proximo[4]))*100);
		#last; ## break o laço for interno
	}
	
	## verifica se o microssatelite atual é maior que o proximo e o overlap é maior que o permitido
	## se sim, não considera a proxima sequencia como microssatelite válido
	if( (($atual[5] - $atual[4]) >= ($proximo[5] - $proximo[4])) && ($overlap_tmp > $overlap) ){
	  
		## elimina o proximo elemento 
		delete $cat_result[$contProximo];
		
		## testa se é o ultimo elemento do hash.
		if($contProximo == $countSequence ){
		      ## guarda a sequencia atual como resultado
		      $result[$qt] = $cat_result[$contAtual];
		      ## incrementa o contador de qtidade de microssatelites válidos
		      $qt++;
		}
	  
		next; ## vai para o proximo elemento do hash
	}
	
	
	## verifica se o microssatelite atual é maior que o proximo e o overlap é maior que o permitido
	## se sim, não considera a proxima sequencia como microssatelite válido
	elsif( (($atual[5] - $atual[4]) < ($proximo[5] - $proximo[4])) && ($overlap_tmp > $overlap) ){
	  
		## elimina o proximo elemento 
		#delete $cat_result[$contProximo];
		
		## testa se é o ultimo elemento do hash.
		#if($contProximo == $countSequence ){
		      ## guarda a sequencia atual como resultado
		#      $result[$qt] = $cat_result[$contAtual];
		      ## incrementa o contador de qtidade de microssatelites válidos
		#      $qt++;
		#}
	  
		#next; ## vai para o proximo elemento do hash
		last; ## incrementa o elemento atual
	}
	##verifica se o microssatelite atual é do mesmo tamanho do proximo 
	##e se o tamanho do overlap está dentro do permitido
	## se sim, salva o microssatelite atual 
	#elsif(($atual[5] - $atual[4]) >= ($proximo[5] - $proximo[4]) && ($overlap_tmp <= $overlap)){
	elsif($overlap_tmp <= $overlap){
		  ## guarda a sequencia atual como resultado
		  $result[$qt] = $cat_result[$contAtual];
		  ## incrementa o contador de qtidade de microssatelites válidos
		  $qt++;
		  ## finaliza o loop interno
		  last;
	}
	##verifica se o microssatelite atual é do mesmo tamanho do proximo 
	##e se o tamanho do overlap é maior que o permitido
	## se sim, excluir o proximo microssatelite
	#elsif(($atual[5] - $atual[4]) < ($proximo[5] - $proximo[4]) && ($overlap_tmp <= $overlap)){
	#	  
	#	  ## elimina o proximo elemento 
	#	 delete $cat_result[$contProximo];
	#	 
	#	 if($contProximo == $countSequence ){
	#	      ## guarda a sequencia atual como resultado
	#	      $result[$qt] = $cat_result[$contAtual];
	#	      ## incrementa o contador de qtidade de microssatelites válidos
	#	      $qt++;
	#	  }
		  
	#}
      }else{
		  ## guarda a sequencia atual como resultado
		  $result[$qt] = $cat_result[$contAtual];
		  ## incrementa o contador de qtidade de microssatelites válidos
		  $qt++;
		  ## finaliza o loop interno
		  last;	  
	}
    } 
  }  
  
  for(my $inc = 0 ; $inc <= $#result; $inc++){
      
      my($pa, $pc, $pg, $pt);
      my $stat = "";
      
      chomp $result[$inc];
      #divide a hash em elementos separados e os armazena em um vetor, onde cada um pode ser acessado separadamente
      my @mostrar = split(/[\t]/, $result[$inc]); 
      my $seq = $db->seq($mostrar[0], $mostrar[4] => $mostrar[5]);
      
      ## Calcula a qtia de cada nucleotide(SOMENTE NUCLEOTIDE)
      if(valida_dna($seq)){
	($pa, $pc, $pg, $pt) = statistic($seq);
	$stat = sprintf ("A(%.2f\%) C(%.2f\%) G(%.2f\%) T(%.2f\%)", $pa,$pc,$pg,$pt);
      }
#      print $mostrar[0]."\t".$mostrar[1]."\t".$mostrar[2]."\t".$mostrar[3]."\t".$mostrar[4]."\t".$mostrar[5]."\t".$mostrar[6]."\t".$mostrar[7]."\n";
      print OUTFILE $mostrar[0]."\t".$mostrar[1]."\t".$mostrar[2]."\t".$mostrar[3]."\t".$mostrar[4]."\t".$mostrar[5]."\t".$mostrar[6]."\t".$mostrar[7]."\t".$stat."\t".$seq."\n";
  }
  $qt = 0;
}

## exclui linha repetidas
system "cat $output_sts | awk '{printf  \"\%s\\n\", \$0}' | sort -u -nk5 > $output_sts.sts";
## deleta arquivos temporarios utilizados
system "rm $output_file";
system "rm $output_file*";
system "rm $output_sts";
##apresenta as informações na tela
print "\nRepeatFinder\t Process has finish. \n";
print "RepeatFinder\t Output filet: $output_sts.sts \n";


sub print_usage {
  
    print <<BLOCK;
USAGE:
   perl progerf.pl -q [FILE] -o [STRING] -i [INTERGER] -y [INTERGER] -r [INTERGER] -g [INTERGER] -v [INTERGER] -d [INTERGER]

DESCRIPTION:
   ProGeRF is a tool for Simple Sequence Repeats (SSR) characterization in any genome or proteome.

REQUIRED ARGUMENTS:
   -q [FILE]
      Input file with DNA/Proteome sequence in FASTA format.
   -o [STRING]
      Output file to be create. default results;
   -i, [INTERGER]
      Minimum length of motif. default 2;
   -y [INTEGER]
      Maximum length of motif. default 5;
   -r [INTERGER] OR -rl [INTERGER-...-INTERGER]
      Minimum repeated times of motif. default 5;   
   -g [INTERGER]
      Maximum allowed Gaps between motifs of a tandem repeat. default 0;
   -v [INTERGER]
      Maximum allowed overllap.  Values allowed  0 - 100. default 0;
   -d [INTERGER]
      Maximum allowed Degeneration motifs. default 0;
   -m [STRING]
      Run mode. n to nucleotideo or p to protein. Default is n.
   

EXAMPLE:
   perl progerf.pl -q Linfantum_JPCM5.fasta -o output -i 2 -y 6 -r 5 -g 3 -v 1 -d 20 -m n

BLOCK
}




## Função que calcular a qtia de cada nucleotides
## em um sequencia
## param <dna_sequence>
## return porcentagem de cada nucletide
sub statistic{
 my $sequence = $_[0];
  
  if(valida_dna($sequence)){
##    print "Sequence valid\n";
  }else{
    print "Sequence invalid\n";
    exit(0);
  }

  my $a=($sequence=~tr/Aa//);
  my $c=($sequence=~tr/Cc//);
  my $g=($sequence=~tr/Gg//);
  my $t=($sequence=~tr/Tt//);
  my $Total=$a+$c+$g+$t;
  
  my $pa =($a/($Total)*100);
  my $pc =($c/($Total)*100);
  my $pg =($g/($Total)*100);
  my $pt =($t/($Total)*100);
  
  return ($pa,$pc,$pg,$pt);
}

##Função que verificar se uma dada sequencia
## contem somente nucleotide
sub valida_dna{
  my $seq = $_[0];
  my $valid = 1;
  if ($seq =~ /[^ACGTacgt]/) {
        my @aux = split (//, $seq);
        foreach my $charac (@aux) {
          if ($charac =~ /[^ACGTacgt\n]/) {
	    $valid = 0;
            last;
          }
        }
   }
   return $valid;
}