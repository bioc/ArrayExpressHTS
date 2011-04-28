#!perl

#
# Andrew Tikhonov  <andrew@ebi.ac.uk>
#

#NT_166433	protein_coding	exon	11955	12166	.	+	.	 gene_id "ENSMUSG00000078424"; transcript_id "ENSMUST00000105219"; exon_number "1"; gene_name "AC007307.2"; transcript_name "AC007307.2-202";
#NT_166433	protein_coding	CDS	12026	12166	.	+	0	 gene_id "ENSMUSG00000078424"; transcript_id "ENSMUST00000105219"; exon_number "1"; gene_name "AC007307.2"; transcript_name "AC007307.2-202"; protein_id "ENSMUSP00000100854";
#NT_166433	protein_coding	start_codon	12026	12028	.	+	0	 gene_id "ENSMUSG00000078424"; transcript_id "ENSMUST00000105219"; exon_number "1"; gene_name "AC007307.2"; transcript_name "AC007307.2-202";
#NT_166433	protein_coding	exon	50322	50492	.	+	.	 gene_id "ENSMUSG00000078424"; transcript_id "ENSMUST00000105219"; exon_number "2"; gene_name "AC007307.2"; transcript_name "AC007307.2-202";
#NT_166433	protein_coding	CDS	50322	50492	.	+	0	 gene_id "ENSMUSG00000078424"; transcript_id "ENSMUST00000105219"; exon_number "2"; gene_name "AC007307.2"; transcript_name "AC007307.2-202"; protein_id "ENSMUSP00000100854";
#NT_166433	protein_coding	exon	51351	51425	.	+	.	 gene_id "ENSMUSG00000078424"; transcript_id "ENSMUST00000105219"; exon_number "3"; gene_name "AC007307.2"; transcript_name "AC007307.2-202";
#NT_166433	protein_coding	CDS	51351	51401	.	+	0	 gene_id "ENSMUSG00000078424"; transcript_id "ENSMUST00000105219"; exon_number "3"; gene_name "AC007307.2"; transcript_name "AC007307.2-202"; protein_id "ENSMUSP00000100854";
#NT_166433	protein_coding	stop_codon	51402	51404	.	+	0	 gene_id "ENSMUSG00000078424"; transcript_id "ENSMUST00000105219"; exon_number "3"; gene_name "AC007307.2"; transcript_name "AC007307.2-202";
#NT_166433	protein_coding	exon	51920	52512	.	+	.	 gene_id "ENSMUSG00000078424"; transcript_id "ENSMUST00000105219"; exon_number "4"; gene_name "AC007307.2"; transcript_name "AC007307.2-202";
#NT_166433	protein_coding	exon	47747	47845	.	+	.	 gene_id "ENSMUSG00000078424"; transcript_id "ENSMUST00000105218"; exon_number "1"; gene_name "AC007307.2"; transcript_name "AC007307.2-201";

use strict;

sub main() {

	if ($#ARGV + 1 < 1) {
		die "\nusage: perl $0 <gtf>\n";
	}

	my $gtffname = $ARGV[0];
	my $outfname = "$gtffname.out";;

	my $total = 0;
	my $loaded = 0;
	my $percentage0 = 0;
	my $percentage1 = 0;

	my $region;
	my $attrs;
	my @zz;
	my $z;

	my @xx;
	my $line; 

	print STDERR "gtf=$gtffname out=$outfname \n";

	print STDERR "processing gtf...\n";

	open(GTF, $gtffname) || die "canot open $gtffname \n";
	open(OUT, ">$outfname") || die "canot open fastqoutfname \n";
	

	$total = -s $gtffname;
	$loaded = 0;

	my $gene_id;
	my $transcript_id;
	my $exon_number;
	my $gene_name;
	my $transcript_name;
	my $protein_id;

	#print OUT "chr\tsource\tfeature\tstart\tend\tscore\tstrand\tphase\tgene_id\ttranscript_id\texon_number\tgene_name\ttranscript_name\tprotein_id\n";

	while($line = <GTF>) {

		$loaded += length($line);

		chomp($line);

		#NT_166433	protein_coding	CDS	51351	51401	.	+	0	 gene_id "ENSMUSG00000078424"; transcript_id "ENSMUST00000105219"; exon_number "3"; gene_name "AC007307.2"; transcript_name "AC007307.2-202"; protein_id "ENSMUSP00000100854";

		#if ($line =~ m/(\S+\t\S+\t\S+\t\S+\t\S+\t\S+\t\S+)\t(.*)^/ ) {
		if ($line =~ m/(.*)\t(.*)$/ ) {
			$region = $1;
			$attrs = $2;

			$gene_id = "";
			$transcript_id = "";
			$exon_number = "";
			$gene_name = "";
			$transcript_name = "";
			$protein_id = "";

			#print STDERR "region=$region\n";
			#print STDERR "attrs=$attrs\n";

			@zz = split(";", $attrs);
			foreach $z (@zz) {
				#print STDERR "orig z=$z\n";
				$z =~ s/ //;
				$z =~ s/\"//g;

				#print STDERR "norm z=$z\n";
	
				@xx = split(" ", $z); 
				
#				gene_id
#				gene_id "ENSMUSG00000078424"; transcript_id "ENSMUST00000105219"; exon_number "3"; gene_name "AC007307.2"; transcript_name "AC007307.2-202"; protein_id "ENSMUSP00000100854";

				#print STDERR "x0=$xx[0]\n";
				#print STDERR "x1=$xx[1]\n";
				
				if ($xx[0] eq "gene_id") { $gene_id = $xx[1]; }
				elsif($xx[0] eq "transcript_id") { $transcript_id = $xx[1]; } 
				elsif($xx[0] eq "exon_number") { $exon_number = $xx[1]; }
				elsif($xx[0] eq "gene_name") { $gene_name = $xx[1]; }
				elsif($xx[0] eq "transcript_name") { $transcript_name = $xx[1]; }
				elsif($xx[0] eq "protein_id") { $protein_id = $xx[1]; } 
				else {
					print STDERR "NO MATCH: $xx[0]\n";
				}
			}

    			#print STDERR "record: gene_id=$gene_id transcript_id=$transcript_id exon_number=$exon_number gene_name=$gene_name transcript_name=$transcript_name protein_id=$protein_id \n";
			print OUT "$region\t$gene_id\t$transcript_id\t$exon_number\t$gene_name\t$transcript_name\t$protein_id\n";

		} else {
			print STDERR "NO MATCH: line=$line\n";
		}

		$percentage0 = int($loaded / $total * 100);

		if ($percentage0 != $percentage1) {
			$percentage1 = $percentage0;
			if ($percentage0 % 10 == 0) {
				print STDERR "processed $percentage0 \%\n";
			}
		}
    	}

	close(OUT);
	close(GTF);

	print STDERR "$outfname saved.\n";
}

main();
1;

