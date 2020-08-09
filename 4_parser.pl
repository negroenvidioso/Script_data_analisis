#!/usr/bin/perl

=for

	Este script parsea los ids de Names_hsa_v22.txt (mirbase) obtenido con 1_GetNames_miRNA_mirbase.pl
	contra las plataformas descargadas (que comienzan con GPL*), 
	de tal forma que se obtienen los microRNA pertenecientes a humano (hsa) que estan siendo estudiados en cada paltaforma.
	
	Ejecutar como:
	ls GPL* | awk 'BEGIN{ORS=" ";}{print $NF;}' | awk '{print "perl 4_parser.pl "$0}' | sh

	Adicionalmente, creara una columna "BINARIA", en el cual si existe el match entre ambos archivos
	se escribira el valor de "1", que significara que "Existe" el microRNA buscada dentro de la plataforma


	Al terminar el script, se crearan los archivos v6_plataforma, donde la columna 3 sera "1 o 0"
	eso indicara si esta presente el microRNA o no, con todos los archivos, crear la matriz binaria:

	ls v6* | awk '{ORS=" ";}{print $NF}' | awk '{print "paste "$0 "| awk \x27{print $1\"\\t\"$3\"\\t\"$6\"\\t\"$9\"\\t\"$12\"\\t\"$15\"\\t\"$18\"\\t\"$21\"\\t\"$24\"\\t\"$27\"\\t\"$30\"\\t\"$33\"\\t\"$36\"\\t\"$39\"\\t\"$42\"\\t\"$45\"\\t\"$48\"\\t\"$51\"\\t\"$54\"\\t\"$57\"\\t\"$60\"\\t\"$63\"\\t\"$66\"\\t\"$69\"\\t\"$72\"\\t\"$75\"\\t\"$78}\x27 > Matriz_Binaria_v6_mirBase22.txt"}' | sh
	
=cut


# Abre directamente el archivo ya que se dara como argumento las plataformas descargadas
my $id_file = "lista_hsa_unique.txt";


# Crea un arreglo con los id del archivo anterior
my $row = `cut -f1 $id_file | awk 'BEGIN{ORS=" ";}{print;}' | less `;
chomp $row;

# Crea un array con los id encontrados en el archivo de mirBase
my @id_hsa = split(" ",$row);

if (scalar(@id_hsa) == 0) {
	system("clear");
	print "\n\t\t*** WARNING ***\n\nNo existe el archivo que contiene los id de los microRNA: $id_file\n";
	print "Buscar el script 1_GetNames_miRNA_mirbase.pl y generarlo en este directorio para ser utilizado\n\n";
	exit;
}
for (my $i = 0 ; $i < scalar(@ARGV) ; $i++ ){
#	Rescata el nombre del archivo leido para usarlo de cabecera despues
	my @nombre_file = split /\./,$ARGV[$i]; 

	my $NombreSalida = "v6_".$ARGV[$i];
	open(FH_1,'>>',$NombreSalida);
	
	print FH_1 "miRNA\t$nombre_file[0]\t$nombre_file[0]\n";
	
	for ( my $j = 0 ; $j < scalar(@id_hsa) ; $j++ ){
		
		my $positivo = 0;
		my $match = 0;

=for ###############################################
		$busqueda
	a) grep -E \"($id_hsa[$j]\\s?)\" $ARGV[$i] 	=> -E = Para buscar expresiones regulares
												=> \"($id_hsa[$j]\\s?)\" = Buscara el ID acompañado de espacio o tabular (en caso de tener, ?: 0 o más)
	b) awk \'END{print NR}\' => Imprimira la cantidad de veces que se encontro el patron anterior
								Ya que grep puede devolver multiples resultados (1 por cada linea encontrada)
=cut ###############################################
		
		my $busqueda = `grep -E \"($id_hsa[$j]\\s?)\" $ARGV[$i] | awk \'END{print NR}\'`;
		
		if ($busqueda){			# Si existe $busqueda
			if($busqueda > 0){ 	# Si $busqueda es mayor que 0
				chomp $busqueda; # Al hacer 2 busquedas, queda con un salto de linea el valor de $busqueda
				$positivo=$busqueda;
				$match = 1;
			}
		}
		else{
			$positivo=0;
		}
		#	Imprimira cuantas veces se encontro el patron buscado en el archivo leido
		print FH_1 "$id_hsa[$j]\t$positivo\t$match\n";
	}
	close(FH_1);
}


