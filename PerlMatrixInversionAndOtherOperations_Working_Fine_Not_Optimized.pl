#! /usr/bin/perl -w

open file1, "<TestDataPerl.csv" || die "Cannot Open File";
binmode file1;

my ($fileContent) = @_;

my $line;

while ($line = readline(file1)) {
	
	next if $line =~ /^\s*$/; 

        if ($line =~ /^([A-Za-z]\w*)/) {

            $fileName = $1;

        } else {
	
	$line =~ s/\r\n$/\n/;
           
	my (@temp) = split(/,/,$line);

	push (@fileMatrix, \@temp);

		}

	}

close file1;

print "\n $fileName \n";

print $fileMatrix[0][0]; 	

my $N = @fileMatrix;

print "\n $N \n";

my $K = @{ $fileMatrix [0] };

print $K;


###############Display Matrix ############################

#display_matrix(@arr);

sub display_matrix (@) {

	my @arr = @_;

	my $N = @arr;

 	my $K = @{ $arr [0] };

	print "\n##################Matrix Display#######################\n";

	@Col = (1 .. $K);

	@row = (1 .. $N);

	$col = join(" ",@Col);

	print "\n\t\t$col\n";

	for($i=0;$i<$N;$i++){

		print "\n$row[$i] ";

		for($j=0;$j<$K;$j++){

			chomp($temp=$arr[$i]->[$j]);

			print " $temp";
			
		}

	}

	}

#display_matrix(@fileMatrix);

sub transpose (@) {

	my @arr = @_;

	my $N = @arr;

 	my $K = @{ $arr [0] };

	for($i=0;$i<$K;$i++){

		for($j=0;$j<$N;$j++){

		$transpose[$i][$j] = $arr[$j][$i];

		}

		}

	return @transpose;			

	}


@transposeMatrix = transpose (@fileMatrix);

display_matrix(@fileMatrix);
print "\n++++++++++++++++++++++++++++++++++\n";

display_matrix(@transposeMatrix);

print "\n+++++++++++++++++++++++++\n";

print @transposeMatrix;

$NK = @transposeMatrix;

$NP = @{ $transposeMatrix[0] };

print "\n N = $NK and K = $NP \n";

##################################Matrix Square###################################


sub matrix_square (@){

	my @arr = @_;

	my $N = @arr;

 	my $K = @{ $arr [0] };

	my @transp = transpose (@arr);

	for($i=0; $i<$K; $i++){

		for($j=0; $j<$K; $j++){
		
			for($p=0; $p<$N; $p++){

			$prod[$i][$j] += ($transp[$i]->[$p])*($arr[$p]->[$j]);

			}
		}
	}

	return @prod;
}


@squareMatrix = matrix_square(@fileMatrix);

print "\n+++++++++++++++++++++++++++++++++\n";

display_matrix(@squareMatrix);

#####################################Matrix Multiply#################################

sub matrix_multiply {

	@arr1 = @{$_[0]};
	@arr2 = @{$_[1]};

	$N1 = @arr1;
	$N2 = @arr2;

	$K1 = @{$arr1[0]};
	$K2 = @{$arr2[0]};

	if($K1 != $N2){
		print "\nERROR: Non conforming matrices. Cannot be multiplied\n";
		return 0;
	}

	for($i=0; $i<$N1; $i++){

		for($j=0; $j<$K2; $j++){

			for($p=0; $p<$K1; $p++){

				$mult[$i][$j] += ($arr1[$i]->[$p])*($arr2[$p]->[$j]);

				}
			}
		
		}

	return @mult;

	}






@multMatrix = matrix_multiply (\@fileMatrix, \@transposeMatrix);

print "\n++++++++++++++++++++++++++++++++++++++++\n";

display_matrix(@multMatrix);


@test1 = (1 .. 7);
@test2 = (3 .. 9);
@test = (\@test1, \@test2);

print "\n++++++++++++++++++++++++++++++++++++++++\n";

matrix_multiply(\@test,\@fileMatrix);

	 	
############################Matrix Addition#################################

sub matrix_addition{

	@arr1 = @{$_[0]};
	@arr2 = @{$_[1]};

	$N1 = @arr1;
	$N2 = @arr2;

	$K1 = @{ $arr1[0] };
	$K2 = @{ $arr2[0] };

	if(($N1 != $N2) || ($K1 != $K2)){
		print "\nERROR: Non conforming matrices. Cannot be added\n";
		return 0;
	}

	for($i=0; $i<$N1; $i++){
		for($j=0; $j<$K1; $j++){
			$add[$i][$j] = ($arr1[$i]->[$j])+($arr2[$i]->[$j]);
		}
	}

	return @add;

	}

print "\n++++++++++++++++++++++++++++++++++++++++\n";


matrix_addition(\@fileMatrix, \@transposeMatrix);

print "\n++++++++++++++++++++++++++++++++++++++++\n";


@addMatrix = matrix_addition(\@fileMatrix, \@fileMatrix);

display_matrix(@addMatrix);

		
############################Matrix Addition#################################

sub matrix_subtraction{

	@arr1 = @{$_[0]};
	@arr2 = @{$_[1]};

	$N1 = @arr1;
	$N2 = @arr2;

	$K1 = @{ $arr1[0] };
	$K2 = @{ $arr2[0] };

	if(($N1 != $N2) || ($K1 != $K2)){
		print "\nERROR: Non conforming matrices. Cannot be added\n";
		return 0;
	}

	for($i=0; $i<$N1; $i++){
		for($j=0; $j<$K1; $j++){
			$subtract[$i][$j] = ($arr1[$i]->[$j])-($arr2[$i]->[$j]);
		}
	}

	return @subtract;

	}

print "\n++++++++++++++++++++++++++++++++++++++++\n";


matrix_subtraction(\@fileMatrix, \@transposeMatrix);

print "\n++++++++++++++++++++++++++++++++++++++++\n";


@subMatrix = matrix_subtraction(\@fileMatrix, \@fileMatrix);

display_matrix(@subMatrix);

##############################Matrix Inversion##################################

sub inverse_matrix (@){

	@arr = @_;

	$N = @arr;

	$K = @{ $arr [0] };

	if($N != $K){

	print "\n ERROR: Not a Square Matrix. Cannot Be inverted \n";

	return 0;

	}

	#Create a conjugate Matrix

	for($i=0; $i<$N; $i++){

		for($j=0; $j<$N; $j++){

			if($i == $j){

			$con[$i][$j] = 1;

			}else{

			$con[$i][$j] = 0;

			}			

		}
	}

	

	for($i=0; $i<($N-1); $i++){
		
		for($j=$i+1; $j<$N; $j++){

			for($k=0; $k<$N; $k++){

				$temp_ji = $arr[$j]->[$i];
				$arr[$j][$k] = ($arr[$j]->[$k])-($temp_ji*($arr[$i]->[$k])/($arr[$i]->[$i]));
				$con[$j][$k] = ($con[$j]->[$k])-($temp_ji*($con[$i]->[$k])/($arr[$i]->[$i]));
					
			}
		
		}
	
	}


	for($i=($N-1); $i>0; $i--){
		
		for($j=$i-1; $j>=0; $j--){

			for($k=0; $k<$N; $k++){

				$temp_ji = $arr[$j]->[$i];
				$arr[$j][$k] = ($arr[$j]->[$k])-($temp_ji*($arr[$i]->[$k])/($arr[$i]->[$i]));
				$con[$j][$k] = ($con[$j]->[$k])-($temp_ji*($con[$i]->[$k])/($arr[$i]->[$i]));
					
			}
		
		}
	
	}

	for($i=0; $i<$N; $i++){

		for($j=0; $j<$N; $j++){

			$tempii = $arr[$i]->[$i];
			$arr[$i][$j] = ($arr[$i]->[$j])/($tempii);
			$con[$i][$j] = ($con[$i]->[$j])/($tempii);			
		
			}

		}	

	display_matrix(@arr);

	return @con;
	
	} 
	

@inverseMatrix = inverse_matrix(matrix_square(@fileMatrix));

display_matrix(@inverseMatrix);

@inverseInverseMatrix = inverse_matrix(@inverseMatrix);

display_matrix(@inverseInverseMatrix);

display_matrix(matrix_square(@fileMatrix));


######################Column Extraction#####################################

sub column_extract(\@$){

	@arr=@{$_[0]};
	$n=$_[1];
	
	$N = @arr;

	$K = @{ $arr [0] };

	if($n > $K){
		print "\nNon COnformable arguments\n";
		return 'NULL';
		}

	for($i=0; $i<$N; $i++){

		chomp($col[$i] = $arr[$i]->[$n]);

		}

	return @col;
		
	

	}

display_matrix(@fileMatrix);

@col=column_extract(@fileMatrix,0);

print "\n#########################################################\n";
print "\n#########################################################\n";
print "\n@col\n";
