#! perl
use strict;

my $count=0;
my $seq;
my @matrix = [];


#$seq = read_data("pattern_matching_data.txt");
#print "Seq: " . substr($seq, 0, 5) . "..." . substr($seq, length($seq) - 5, 5) . "\n";

###pattern count
#$count = pattern_count($seq,"ACACCA");
#print "Pattern appears " . $count . " times\n";

###find all the k-mers in a sequence
#read_kmers($seq, 12);

###get dna reverse complement 
#traverse_reverse_seq($seq);

###calculate and show the skew of a sequence
#show_skew("TAAAGACTGCCGAGAGGCCAACACGAGTGCTAGAACGAGGGGCGTAAACGCGGGTCCGAT");

###calculate hamming distance
#my $hamming_dist = calculate_hamming_distance("GGGCCGTTGGT", "GGACCGTTGAC");
#print $hamming_dist;

hamming_sequence("CGCCCGAATCCAGAACGCATTCCCATATTTCGGGACCACTGGCCTCCACGGTACGGACGTCAATCAAAT", "ATTCTGGA", 3);

sub hamming_sequence()
{
	my $text = uc($_[0]);
	my $pattern = uc($_[1]);
	my $max_hdist = uc($_[2]);
	my $i;

	for($i=0; $i <= length($text) - length($pattern); $i++)
	{
		if (calculate_hamming_distance(substr($text, $i, length($pattern)), $pattern) <= $max_hdist)
		{
			print "$i ";
		}
	}
}

sub calculate_hamming_distance()
{
	my $text = $_[0];
	my $pattern = $_[1];

	my $hdist = 0;

	for(my $i=0; $i<length($text); $i++)
	{
		if (index(substr($text,$i,1), substr($pattern,$i,1)) < 0)
		{
			$hdist++;
		}
	}
	return $hdist;
}

sub show_skew()
{
	my $text = $_[0];
	my $count_G = 0;
	my $count_C = 0;
	my $delta = 0;
	my @skew_array;
	my $skew_index = 0;
	my $min;

	print $delta . " ";
	$skew_array[$skew_index] = 0;
	$skew_index++;
	while($text =~ /(.)/g)
	{
		if (index($1, "G") == 0)
		{
			$count_G++
		}
		if (index($1, "C") == 0)
		{
			$count_C++;
		}
		$delta = $count_G - $count_C;
		$skew_array[$skew_index] = $delta;
		$skew_index++;
		print $delta . " ";
	}
	$min = find_min_skew(\@skew_array, $skew_index);
	print "Min: " . $min . "\n";
	show_min_skew(\@skew_array, $skew_index, $min);
}

sub find_min_skew()
{
	my @skew_array = @{$_[0]};
	my $n = $_[1];
	my $min = 0;

	print "n=" . $n . "\n";
	for(my $i=0; $i<$n;$i++)
	{
		#print "check " . $skew_array[$i] . "\n";
		if ($skew_array[$i] < $min)
		{
			$min = $skew_array[$i];
		}
	}
	return $min;
}

sub show_min_skew()
{
	my @skew_array = @{$_[0]};
	my $n = $_[1];
	my $min = $_[2];

	for(my $i=0; $i<$n;$i++)
	{
		if ($skew_array[$i] == $min)
		{
			print $i . " ";
		}
	}
}

sub traverse_reverse_seq()
{
	my $text = $_[0];
	my $dna_rev = "";

	for(my $i=length($text)-1; $i>=0; $i--)
	{
		$dna_rev .= dna_complement(substr($text, $i, 1));
	}
	print "DNA complement: \n" . $dna_rev . "\n";
}

#params: nucleotide(amino acid)
sub dna_complement()
{
	my $aa = uc($_[0]);
	if ($aa =~ /A/)
	{
		return "T";
	}
	if ($aa =~ /T/)
	{
		return "A";
	}
	if ($aa =~ /C/)
	{
		return "G";
	}	
	if ($aa =~ /G/)
	{
		return "C";
	}
}

#params: text, kmer_len
sub read_kmers()
{
	my $text = uc($_[0]);
	my $kmer_len = uc($_[1]);
	my $i;
	my $pattern = "";
	my $index = 0;
	my $max;

	for($i=0; $i <= length($text) - $kmer_len; $i++)
	{
		$pattern = substr($text, $i, $kmer_len);
		#print "Get pattern: " . $pattern . "\n";
		if (is_pattern_in_matrix($index, $pattern) == 0)
		{
			add_pattern_to_matrix($index, $pattern, pattern_count($text, $pattern));
			$index++;
		}
	}
	$max = get_max_kmer($index);
	print "Max kmer: " . $max . "\n";
	show_max_kmers($index, $max);
}

sub is_pattern_in_matrix()
{
	my $index = $_[0];
	my $pattern = $_[1];

	for(my $i=0; $i<$index; $i++)
	{
		if (index($matrix[$i][0], $pattern) == 0)
		{
			return 1;
		}
	}
	return 0;
}

#params: matrix, index, pattern, pattern count
sub add_pattern_to_matrix()
{
	my $index 	= $_[0];
	my $pattern = $_[1];
	my $pcount 	= $_[2];

	$matrix[$index][0] = $pattern;
	$matrix[$index][1] = $pcount;
	#print "Added " . $matrix[$index][0] . " with count " . $matrix[$index][1] . "\n";
}

sub get_max_kmer()
{
	my $index	= $_[0];
	my $max_kmer=0;

	for(my $i=0; $i<$index; $i++)
	{
		#print "Check " . $matrix[$i][0] . "\n";
		if ($matrix[$i][1] > $max_kmer)
		{
			$max_kmer = $matrix[$i][1];
		}
	}
	return $max_kmer;
}

sub show_max_kmers()
{
	my $index	= $_[0];
	my $max_kmer= $_[1];
	for(my $i=0; $i<$index; $i++)
	{
		if ($matrix[$i][1] == $max_kmer)
		{
			print "\t" . $matrix[$i][0];
		}
	}
}

#params: text, pattern
sub pattern_count()
{
	my $text = uc($_[0]);
	my $pattern = uc($_[1]);
	my $i;
	my $tsize = length($text);
	my $psize = length($pattern);
	my $aux;
	my $counter = 0;

	if ($tsize == 0 || $psize == 0 || $tsize < $psize)
	{
		return $counter;
	}

	for($i=0; $i<=$tsize - $psize; $i++)
	{
		$aux = substr($text,$i, $psize);
		if (index($aux, $pattern) == 0)
		{
			print "found @ " . $i . "\n";
			$counter++;
		}
	}
	return $counter;
}

sub read_data()
{
	my $filename = $_[0];
	my $text = "";
	my $flag = 1;
	open(IN, "<", $filename);
	while(<IN>)
	{
		chomp();
		if ($_ =~ /Output/)
		{
			$flag = 0;
		}
		if ($_ =~ /([ATCG]+)/)
		{
			if (length($1) > 50 && $flag == 1)
			{
				$text .= uc($1);
			}
		}
	}
	close(IN);
	return $text;
}