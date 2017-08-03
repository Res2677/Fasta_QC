use strict;
use FindBin '$Bin';
use Getopt::Long;
use List::Util qw/max min/;
use Thread;
use Test::More;
use Data::Dumper;
my ($file1,$file2,$outdir,$help,$seg,$Rscript,$Convert,$seg_4);
GetOptions(
        "f1=s" => \$file1,
        "f2:s" => \$file2,
        "o=s"  => \$outdir,
	"seg:i" => \$seg,
        "help|?" => \$help
);

if (!defined $file1 || !defined $outdir ||defined $help) {
        die << "USAGE";
        description: cluster
        usage: perl $0 [options]
options:
        -f1 <file>    *fq.gz (SE,or PE:read1)
        -f2 <file>    *fq.gz (PE read2)
        -o  <path>    path of outdir
	-seg  number   number of files that have been split up(The more the number, the faster the speed,default: 20)
        -help|?       information for help
e.g.:
perl $0 -f1 /share/stor-01/zebra-AUO/01.WGS/liangxinming/20160914.NA12878_PE100/01.filter/CL100008328_L01/test_1.fq.gz -f2 /share/stor-01/zebra-AUO/01.WGS/liangxinming/20160914.NA12878_PE100/01.filter/CL100008328_L01/test_2.fq.gz -o /ifs4/BC_RD/USER/huangjingfeng/2017-2-06/QC_BOX/TEST
USAGE
}
$Rscript    ||= "/ifs4/BC_PUB/biosoft/pipeline/Package/R-3.1.1/lib64/R/bin/Rscript";
$Convert    ||= "/usr/bin/convert";
$seg	    ||= 10;

$seg_4 = $seg*4;
system("echo ready: `date`");
system("mkdir -p $outdir/process; mkdir -p $outdir/result");
my($decom_fqgz,$decom_fq);
my($seqx4_1,$seqx4_2);
my($split1,$split2,$seq_num_1,$seq_num_2);
if (!defined $file2)
{
	if ($file1=~/.fq.gz$/) {$seqx4_1 = &decom_fqgz($file1,"read1");}
	elsif ($file1=~/.fq$/) {$seqx4_1 = &decom_fq($file1,"read1");}
	else {print "error : cannot identify this datatype";}
}
else
{
        my($split1,$split2);
        if ($file1=~/.fq.gz$/) {$split1 = Thread->new($seqx4_1 = &decom_fqgz($file1,"read1"));}
	elsif ($file1=~/.fq$/) {$split1 = Thread->new($seqx4_1 = &decom_fq($file1,"read1"));}
	else {print "read1 error : cannot identify this datatype";}
        if ($file2=~/.fq.gz$/) {$split2 = Thread->new($seqx4_1 = &decom_fqgz($file2,"read2"));}
	elsif ($file2=~/.fq$/) {$split2 = Thread->new($seqx4_1 = &decom_fq($file2,"read2"));}
	else {print "read2 error : cannot identify this datatype";}
        $split1 -> join();
        $split2 -> join();
}
print "$seqx4_1\t$seqx4_2\n";
$seq_num_1 = $seqx4_1/4;
if (defined $file2){$seq_num_2 = $seqx4_2/4;}

sub decom_fqgz
{
        my $file = shift;
        my $read = shift;
        my $word = `gzip -dc $file > $outdir/result/$read.fq; echo gzip: \`date\`; WORD=\`cat $outdir/result/$read.fq|wc -l|awk '{print int(\$1/$seg_4)*4,\$1}'\`;read -a ARR <<< \$WORD;echo \${ARR[1]}; split -l \${ARR[0]} $outdir/result/$read.fq $outdir/process/$read.fq; echo split: \`date\``;
        my @ww = split /\n/,$word;
        return $ww[1];
}

sub decom_fq
{
	my $file = shift;
	my $read = shift;
	my $word = `WORD=\`cat $file |wc -l|awk '{print int(\$1/$seg_4)*4,\$1}'\`;read -a ARR <<< \$WORD;echo \${ARR[1]}; split -l \${ARR[0]} $file $outdir/process/$read.fq; echo split: \`date\``;
	my @ww = split /\n/,$word;
	return $ww[1];
}

my (@name1,$file1_name,@per_qc,@name2,$file2_name,@a,@file_seg1,$file1_num,@file_seg2,$file2_num,$per_read_num,$base_line,$qual_line,$read_length,$q,$min_qc,$Encoding,$asi_qc,$nn,@ss,@ss1,@ss2,$pos,$v,$qual,$a_content,$t_content,$c_content,$g_content,$n_content,$seq_gc,$content,$number,%sequence_gc,%qual_count_pos,%qual_count,$Q0,$qc,$q_num,$lower_qual,$smaller_qual_ranking,$smaller_qual_ranking,$t,$c,$g,$n,$mean);
@file_seg1 = glob ("$outdir/process/read1.fq*");$file1_num = @file_seg1;
if (defined $file2)
{
	@file_seg2 = glob ("$outdir/process/read2.fq*");$file2_num = @file_seg2;
}
#print Dumper \@file_seg2;
@name1 = split/\//,$file1;
$file1_name = $name1[-1];
if (defined $file2)
{
	@name2 = split/\//,$file2;
	$file2_name = $name2[-1];
}

open FQ1,"$file_seg1[0]",or die;
while(<FQ1>)
{
    chomp;
    $per_read_num ++;
    if ($per_read_num == 100){last;}
    $base_line = <FQ1>;
    <FQ1>;
    $qual_line = <FQ1>;
    chomp $base_line;
    chomp $qual_line;
    my @qual = split //,$qual_line;
    my @base = split //,$base_line;
    $read_length = @qual;
    foreach $q (@qual)
    {
        $q = ord $q;
        push @per_qc,$q;
    }
}
close FQ1;
$min_qc = min @per_qc;
if ($min_qc < 64){$Encoding = "phred33";$asi_qc = 33;}
elsif ($min_qc >= 64){$Encoding = "phred64";$asi_qc = 64;}
my $seq_num = $seq_num_1 + $seq_num_2;
our $date_cnt = ($seq_num*$read_length)/1000000000;
system("echo \"start (base_num: $date_cnt G): `date`\"");

open POUT1,">$outdir/process/read1.Basic_Statistics.xls",or die $!;
print POUT1 "GC%\tN%\tQ20\tQ30\n";
close POUT1;
open POUT2,">$outdir/process/read1.Per_base_sequence_quality.xls",or die $!;
print POUT2 "Position\ttype\value\n";
close POUT2;
open POUT3,">$outdir/process/read1.Per_sequence_quality_scores.xls",or die $!;
print POUT3 "qc\tnumber\n";
close POUT3;
open POUT4,">$outdir/process/read1.Per_base_sequence_content.xls",or die $!;
print POUT4 "pos\tbase\tnumber\n";
close POUT4;
open POUT5,">$outdir/process/read1.Per_sequence_GC_content.xls",or die $!;
print POUT5 "content\tnumber\n";
close POUT5;
open POUT1,">>$outdir/process/read1.Basic_Statistics.xls",or die $!;
open POUT2,">>$outdir/process/read1.Per_base_sequence_quality.xls",or die $!;
open POUT3,">>$outdir/process/read1.Per_sequence_quality_scores.xls",or die $!;
open POUT4,">>$outdir/process/read1.Per_base_sequence_content.xls",or die $!;
open POUT5,">>$outdir/process/read1.Per_sequence_GC_content.xls",or die $!;
if(defined $file2)
{
open ROUT1,">$outdir/process/read2.Basic_Statistics.xls",or die $!;
print ROUT1 "seq_num\tGC%\tN%\tQ20\tQ30\n";
close ROUT1;
open ROUT2,">$outdir/process/read2.Per_base_sequence_quality.xls",or die $!;
print ROUT2 "Position\ttype\value\n";
close ROUT2;
open ROUT3,">$outdir/process/read2.Per_sequence_quality_scores.xls",or die $!;
print ROUT3 "qc\tnumber\n";
close ROUT3;
open ROUT4,">$outdir/process/read2.Per_base_sequence_content.xls",or die $!;
print ROUT4 "pos\tbase\tnumber\n";
close ROUT4;
open ROUT5,">$outdir/process/read2.Per_sequence_GC_content.xls",or die $!;
print ROUT5 "content\tnumber\n";
close ROUT5;
open ROUT1,">>$outdir/process/read2.Basic_Statistics.xls",or die $!;
open ROUT2,">>$outdir/process/read2.Per_base_sequence_quality.xls",or die $!;
open ROUT3,">>$outdir/process/read2.Per_sequence_quality_scores.xls",or die $!;
open ROUT4,">>$outdir/process/read2.Per_base_sequence_content.xls",or die $!;
open ROUT5,">>$outdir/process/read2.Per_sequence_GC_content.xls",or die $!;
}
#our ($seq_num1,$GC1,$N_num1,%base1,%g_num1,%c_num1,%a_num1,%t_num1,%n_num1,%sequence_gc1,%qual_count1);
#my %b	: shared;my %c   : shared;
#our %qual_count_pos1;
#our $seq_num1 = 1;

#my($seq_num2,$GC2,$N_num2,%base2,%g_num2,%c_num2,%a_num2,%t_num2,%n_num2,%sequence_gc2,%qual_count_pos2,%qual_count2)   : shared;
#
if (!defined $file2)
{
	foreach $nn (0..$file1_num-1)
	{
		$ss[$nn] = Thread->new(\&start_thread,$file_seg1[$nn],"read1",$seq_num_1);
	}
	foreach $nn (0..$file1_num-1)
	{
        	$ss[$nn] -> join();
	}
	
}
else
{
	foreach $nn (0..$file1_num-1)
	{
		$ss1[$nn] = Thread->new(\&start_thread,$file_seg1[$nn],"read1",$seq_num_1);
	}
	foreach $nn (0..$file2_num-1)
	{
		$ss2[$nn] = Thread->new(\&start_thread,$file_seg2[$nn],"read2",$seq_num_2);
	}
	foreach $nn (0..$file1_num-1)
	{
        	$ss1[$nn] -> join();
	}
	foreach $nn (0..$file2_num-1)
	{
		$ss2[$nn] -> join();
	}
}
close POUT1;
close POUT2;
close POUT3;
close POUT4;
close POUT5;
if (defined $file2)
{
close ROUT1;
close ROUT2;
close ROUT3;
close ROUT4;
close ROUT5;
}

system("echo end: `date`");

my ($GC_per,$N_per,$Q20,$Q30,%Mean,%Median,%Lower_Quartile,%Upper_Quartile,%Percentile_10th,%Percentile_90th,%qual_num,%base_content,%gc_content);
open IN1,"$outdir/process/read1.Basic_Statistics.xls",or die $!;
<IN1>;
while(<IN1>)
{
	chomp;
	@a = split /\t/;
	$GC_per = $GC_per + $a[0];
	$N_per = $N_per + $a[1];
	$Q20 = $Q20 + $a[2];
	$Q30 = $Q30 + $a[3];	
}
$GC_per = sprintf ("%.2f",$GC_per);
$N_per = sprintf ("%.2f",$N_per);
$Q20 = sprintf ("%.2f",$Q20);
$Q30 = sprintf ("%.2f",$Q30);
#$base_num = $seq_num * $read_length;
close IN1;
open OUT1,">$outdir/result/read1.Basic_Statistics.xls";
print OUT1 "Measure\tValue\n";
close OUT1;
open OUT1,">>$outdir/result/read1.Basic_Statistics.xls";
print OUT1 "Filename\t$file1_name\n";
print OUT1 "Encoding\t$Encoding\n";
print OUT1 "Total Sequences\t$seq_num_1\n";
print OUT1 "Sequence length\t$read_length\n";
print OUT1 "%GC\t$GC_per%\n%N\t$N_per%\n";
print OUT1 "%Q20\t$Q20%\n%Q30\t$Q30%\n";
close OUT1;

open IN2,"$outdir/process/read1.Per_base_sequence_quality.xls";
<IN2>;
while(<IN2>)
{
	chomp;
	@a = split /\t/;
	if($a[1] eq "Mean"){$Mean{$a[0]} = $Mean{$a[0]} + $a[2];}
	elsif($a[1] eq "Median"){$Median{$a[0]}{$a[2]} ++;}
	elsif($a[1] eq "Lower Quartile"){$Lower_Quartile{$a[0]}{$a[2]} ++;}
	elsif($a[1] eq "Upper Quartile"){$Upper_Quartile{$a[0]}{$a[2]} ++;}
	elsif($a[1] eq "10th Percentile"){$Percentile_10th{$a[0]}{$a[2]} ++;}
	elsif($a[1] eq "90th Percentile"){$Percentile_90th{$a[0]}{$a[2]} ++;}	
}
close IN2;
open OUT2,">$outdir/result/read1.Per_base_sequence_quality.xls";
print OUT2 "Position\tMean\tMedian\tLower Quartile\tUpper Quartile\t10th Percentile\t90th Percentile\n";
close OUT2;
open OUT2,">>$outdir/result/read1.Per_base_sequence_quality.xls";
foreach $pos (sort {$a <=>$b} keys %Median)
{
	my ($median,$lower_Quartile,$upper_Quartile,$percentile_10th,$percentile_90th);
	foreach $qual (sort {$Median{$pos}{$b} <=> $Median{$pos}{$a}} (keys %{$Median{$pos}})) {$median = $qual;last;}
	foreach $qual (sort {$Lower_Quartile{$pos}{$b} <=> $Lower_Quartile{$pos}{$a}} (keys %{$Lower_Quartile{$pos}})){$lower_Quartile = $qual;last;}
	foreach $qual (sort {$Upper_Quartile{$pos}{$b} <=> $Upper_Quartile{$pos}{$a}} (keys %{$Upper_Quartile{$pos}})){$upper_Quartile = $qual;last;}
	foreach $qual (sort {$Percentile_10th{$pos}{$b} <=> $Percentile_10th{$pos}{$a}} (keys %{$Percentile_10th{$pos}})){$percentile_10th = $qual;last;}
	foreach $qual (sort {$Percentile_90th{$pos}{$b} <=> $Percentile_90th{$pos}{$a}} (keys %{$Percentile_90th{$pos}})){$percentile_90th = $qual;last;}
	$mean = $Mean{$pos}/$file1_num;
	print OUT2 "$pos\t$mean\t$median\t$lower_Quartile\t$upper_Quartile\t$percentile_10th\t$percentile_90th\n";
}
close OUT2;
open IN3,"$outdir/process/read1.Per_sequence_quality_scores.xls";
<IN3>;
while(<IN3>)
{
	chomp;
        @a = split /\t/;
	$qual_num{$a[0]} = $qual_num{$a[0]} + $a[1];
}
close IN3;
open OUT3,">$outdir/result/read1.Per_sequence_quality_scores.xls";
print OUT3 "qc\tnumber\n";
close OUT3;
open OUT3,">>$outdir/result/read1.Per_sequence_quality_scores.xls";
foreach $qual (sort {$a <=> $b} keys %qual_num)
{
	print OUT3 "$qual\t$qual_num{$qual}\n";
}
close OUT3;
open IN4,"$outdir/process/read1.Per_base_sequence_content.xls";
<IN4>;
while(<IN4>)
{
	chomp;
        @a = split /\t/;
	$base_content{$a[0]}{$a[1]} = $base_content{$a[0]}{$a[1]} + $a[2];
}
close IN4;
open OUT4,">$outdir/result/read1.Per_base_sequence_content.xls";
print OUT4 "pos\tA\tT\tC\tG\tN\n";
close OUT4;
open OUT4,">>$outdir/result/read1.Per_base_sequence_content.xls";
foreach $pos (sort {$a <=> $b} keys %base_content)
{
	$a_content = $base_content{$pos}{"A"};
	$t_content = $base_content{$pos}{"T"};
	$c_content = $base_content{$pos}{"C"};
	$g_content = $base_content{$pos}{"G"};
	$n_content = $base_content{$pos}{"N"};
	print OUT4 "$pos\t$a_content\t$t_content\t$c_content\t$g_content\t$n_content\n";
}
close OUT4;
open IN5,"$outdir/process/read1.Per_sequence_GC_content.xls";
<IN5>;
while(<IN5>)
{
	chomp;
        @a = split /\t/;
	$gc_content{$a[0]} = $gc_content{$a[0]} + $a[1];
}
close IN5;
open OUT5,">$outdir/result/read1.Per_sequence_GC_content.xls";
print OUT5 "gc\tnumber\n";
close OUT5;
open OUT5,">>$outdir/result/read1.Per_sequence_GC_content.xls";
foreach $seq_gc (sort {$a <=> $b} keys %gc_content)
{
	$content = ($seq_gc/$read_length)*100;
	$number = $gc_content{$seq_gc};
	print OUT5 "$content\t$number\n";
}
close OUT5;

#############for read2
if (defined $file2)
{
my($GC_per,$N_per,$Q20,$Q30,%Mean,%Median,%Lower_Quartile,%Upper_Quartile,%Percentile_10th,%Percentile_90th,%qual_num,%base_content,%gc_content);
open IN1,"$outdir/process/read2.Basic_Statistics.xls",or die $!;
<IN1>;
while(<IN1>)
{
	chomp;
	@a = split /\t/;
	$GC_per = $GC_per + $a[0];
	$N_per = $N_per + $a[1];
	$Q20 = $Q20 + $a[2];
	$Q30 = $Q30 + $a[3];	
}
$GC_per = sprintf ("%.2f",$GC_per);
$N_per = sprintf ("%.2f",$N_per);
$Q20 = sprintf ("%.2f",$Q20);
$Q30 = sprintf ("%.2f",$Q30);
close IN1;
open OUT1,">$outdir/result/read2.Basic_Statistics.xls";
print OUT1 "Measure\tValue\n";
close OUT1;
open OUT1,">>$outdir/result/read2.Basic_Statistics.xls";
print OUT1 "Filename\t$file2_name\n";
print OUT1 "Encoding\t$Encoding\n";
print OUT1 "Total Sequences\t$seq_num_2\n";
print OUT1 "Sequence length\t$read_length\n";
print OUT1 "%GC\t$GC_per%\n%N\t$N_per%\n";
print OUT1 "%Q20\t$Q20%\n%Q30\t$Q30%\n";
close OUT1;

open IN2,"$outdir/process/read2.Per_base_sequence_quality.xls";
<IN2>;
while(<IN2>)
{
	chomp;
	@a = split /\t/;
	if($a[1] eq "Mean"){$Mean{$a[0]} = $Mean{$a[0]} + $a[2];}
	elsif($a[1] eq "Median"){$Median{$a[0]}{$a[2]} ++;}
	elsif($a[1] eq "Lower Quartile"){$Lower_Quartile{$a[0]}{$a[2]} ++;}
	elsif($a[1] eq "Upper Quartile"){$Upper_Quartile{$a[0]}{$a[2]} ++;}
	elsif($a[1] eq "10th Percentile"){$Percentile_10th{$a[0]}{$a[2]} ++;}
	elsif($a[1] eq "90th Percentile"){$Percentile_90th{$a[0]}{$a[2]} ++;}	
}
close IN2;
open OUT2,">$outdir/result/read2.Per_base_sequence_quality.xls";
print OUT2 "Position\tMean\tMedian\tLower Quartile\tUpper Quartile\t10th Percentile\t90th Percentile\n";
close OUT2;
open OUT2,">>$outdir/result/read2.Per_base_sequence_quality.xls";
foreach $pos (sort {$a <=> $b} keys %Median)
{
	my ($mean,$median,$lower_Quartile,$upper_Quartile,$percentile_10th,$percentile_90th);
	foreach $qual (sort {$Median{$pos}{$b} <=> $Median{$pos}{$a}} (keys %{$Median{$pos}})) {$median = $qual;last;}
	foreach $qual (sort {$Lower_Quartile{$pos}{$b} <=> $Lower_Quartile{$pos}{$a}} (keys %{$Lower_Quartile{$pos}})){$lower_Quartile = $qual;last;}
	foreach $qual (sort {$Upper_Quartile{$pos}{$b} <=> $Upper_Quartile{$pos}{$a}} (keys %{$Upper_Quartile{$pos}})){$upper_Quartile = $qual;last;}
	foreach $qual (sort {$Percentile_10th{$pos}{$b} <=> $Percentile_10th{$pos}{$a}} (keys %{$Percentile_10th{$pos}})){$percentile_10th = $qual;last;}
	foreach $qual (sort {$Percentile_90th{$pos}{$b} <=> $Percentile_90th{$pos}{$a}} (keys %{$Percentile_90th{$pos}})){$percentile_90th = $qual;last;}
	$mean = $Mean{$pos}/$file1_num;
	print OUT2 "$pos\t$mean\t$median\t$lower_Quartile\t$upper_Quartile\t$percentile_10th\t$percentile_90th\n";
}
close OUT2;
open IN3,"$outdir/process/read2.Per_sequence_quality_scores.xls";
<IN3>;
while(<IN3>)
{
	chomp;
        @a = split /\t/;
	$qual_num{$a[0]} = $qual_num{$a[0]} + $a[1];
}
close IN3;
open OUT3,">$outdir/result/read2.Per_sequence_quality_scores.xls";
print OUT3 "qc\tnumber\n";
close OUT3;
open OUT3,">>$outdir/result/read2.Per_sequence_quality_scores.xls";
foreach $qual (sort {$a <=>$b} keys %qual_num)
{
	print OUT3 "$qual\t$qual_num{$qual}\n";
}
close OUT3;
open IN4,"$outdir/process/read2.Per_base_sequence_content.xls";
<IN4>;
while(<IN4>)
{
	chomp;
        @a = split /\t/;
	$base_content{$a[0]}{$a[1]} = $base_content{$a[0]}{$a[1]} + $a[2];
}
close IN4;
open OUT4,">$outdir/result/read2.Per_base_sequence_content.xls";
print OUT4 "pos\tA\tT\tC\tG\tN\n";
close OUT4;
open OUT4,">>$outdir/result/read2.Per_base_sequence_content.xls";
foreach $pos (sort {$a <=> $b} keys %base_content)
{
	$a_content = $base_content{$pos}{"A"};
	$t_content = $base_content{$pos}{"T"};
	$c_content = $base_content{$pos}{"C"};
	$g_content = $base_content{$pos}{"G"};
	$n_content = $base_content{$pos}{"N"};
	print OUT4 "$pos\t$a_content\t$t_content\t$c_content\t$g_content\t$n_content\n";
}
close OUT4;
open IN5,"$outdir/process/read2.Per_sequence_GC_content.xls";
<IN5>;
while(<IN5>)
{
	chomp;
        @a = split /\t/;
	$gc_content{$a[0]} = $gc_content{$a[0]} + $a[1];
}
close IN5;
open OUT5,">$outdir/result/read2.Per_sequence_GC_content.xls";
print OUT5 "gc\tnumber\n";
close OUT5;
open OUT5,">>$outdir/result/read2.Per_sequence_GC_content.xls";
foreach $seq_gc (sort {$a <=> $b} keys %gc_content)
{
	$content = ($seq_gc/$read_length)*100;
	$number = $gc_content{$seq_gc};
	print OUT5 "$content\t$number\n";
}
close OUT5;
}

if (!defined $file2)
{
open my $fh_rcode,">","$outdir/result/SE.QC_Statistics.R" or die $!;
print $fh_rcode <<CMD;
pdf("$outdir/result/SE.QC_Statistics.pdf",8,5)
par(mfrow = c(1,1))
data1 <- read.table("$outdir/result/read1.Basic_Statistics.xls",head=T,sep="\\t")
len <- length(row.names(data1))
xc <- seq(0,len+2,1)
yc <- seq(0,len+2,1)
plot(xc,yc,ylab="",xlab="",cex.lab=1.5,yaxt="n",xaxt="n",axes=F,main="",cex=0)
cc1 <- seq(0.5,len+0.5,0.05)
cc2 <- seq(len+0.5,len+1.5,0.05)
for (i in cc1)
{
        segments(0,i,len+2,i,lwd=4,col="#DEEBF7")
}
for (i in cc2)
{
        segments(0,i,len+2,i,lwd=4,col="#9ECAE1")
}
segments(0,len+0.5,len+2,len+0.5,lwd=3,col="grey20")
segments(0,len+1.5,len+2,len+1.5,lwd=3,col="grey20")
segments(0,0.5,len+2,0.5,lwd=3,col="grey20")
mtext(side=3,"Basic Statistics",line=0,cex=1.8)
text(3,length(row.names(data1))+1,labels="Measure")
text(7,length(row.names(data1))+1,labels="Value")
for (i in 1:length(row.names(data1)))
{
        text(3,(length(row.names(data1))+1- i),labels=data1\$Measure[i])
        text(7,(length(row.names(data1))+1- i),labels=data1\$Value[i])
}
#########1.Basic_Statistics :end

par(mfrow = c(1,1))
data2 <- read.table("$outdir/result/read1.Per_base_sequence_quality.xls",head=T,sep="\\t")
pos_len <- length(row.names(data2))
k <- seq(0,40-40/pos_len,40/pos_len)
w0 <- seq(0,20,0.05)
w20 <- seq(20,30,0.05)
w30 <- seq(30,40,0.05)
par(mar=c(4,5,4,2),xpd=TRUE)
plot(data2\$Position,k,xlab="",ylab="",xaxt="n",yaxt="n",cex=0,ylim = c(0,40),axes = FALSE)
for (i in 1:length(w0))
{
        segments(0,w0[i],max(data2\$Position)+1,w0[i],lwd=2,col=adjustcolor("#C6DBEF",alpha.f = 0.2))
}
for (i in 1:length(w20))
{
        segments(0,w20[i],max(data2\$Position)+1,w20[i],lwd=2,col=adjustcolor("#6BAED6",alpha.f = 0.2))
}
for (i in 1:length(w30))
{
        segments(0,w30[i],max(data2\$Position)+1,w30[i],lwd=2,col=adjustcolor("#2171B5",alpha.f = 0.2))
}
par(new=T)
par(mar=c(4,5,4,2),xpd=TRUE)
xlabel <-c(1,seq(10,max(data2\$Position),10))
ylabel <- seq(0,40,10)
seg.back <- seq(-0.3,0.3,0.02)
plot(data2\$Position,k,xlab="",ylab="",xaxt="n",yaxt="n",cex=0,ylim = c(0,40),axes = FALSE)
for (i in 1:pos_len)
{
        segments(i-0.2,data2\$X10th.Percentile[i],i+0.2,data2\$X10th.Percentile[i],lwd=1,col="grey50")
        segments(i-0.2,data2\$X90th.Percentile[i],i+0.2,data2\$X90th.Percentile[i],lwd=1,col="grey50")
        segments(i,data2\$X10th.Percentile[i],i,data2\$X90th.Percentile[i],lwd=1,col="grey50",lty=2)
	for (j in seg.back)
        {
                segments(i+j,data2\$Lower.Quartile[i],i+j,data2\$Upper.Quartile[i],lwd=0.5,col="#FB9A99")
        }
        segments(i-0.3,data2\$Median[i],i+0.3,data2\$Median[i],lwd=1,col=adjustcolor("#E41A1C",alpha.f = 0.8))
}
lines(data2\$Mean, lty =1, col = "#6A51A3")
axis(1,lwd=2,col="grey40",col.axis="grey40",pos=0,at=xlabel,labels =xlabel)
axis(2,lwd=2,col="grey40",col.axis="grey40",pos=0,at=ylabel,label=ylabel,las=2)
segments(0,0,0,40,lwd=2)
segments(0,0,max(data2\$Position)+1,0,lwd=2)
segments(0,20,max(data2\$Position)+1,20,lwd=2,col="grey70",lty=2)
segments(0,30,max(data2\$Position)+1,30,lwd=2,col="grey70",lty=2)
mtext(side=2,"Quality Scores",line=3,cex=1.5)
mtext(side=3,"Per base sequence quality(read1)",line=1,cex=1.8)
mtext(side=1,"Position in read(bp)",line=2,cex=1.5)
######################2.Per_base_sequence_quality: end

par(mfrow = c(1,1))
data3 <- read.table("$outdir/result/read1.Per_sequence_quality_scores.xls",head=T,sep = "\\t")
xlabel <- c(min(data3\$qc),seq(10,40,5))
par(mar=c(4,5,6,4),xpd=TRUE)
plot(data3\$qc,data3\$number,ylab="",xlab="",xaxt="n",ylim=c(0,max(data3\$number)*1.2),xlim=c(min(data3\$qc),38),axes=F,cex=0)
lines(data3\$qc,data3\$number,lty=1,col="#dc143c",lwd=2)
par(new=T)
par(mar=c(4,5,6,4),xpd=TRUE)
segments(min(data3\$qc),0,max(data3\$qc)+1,0,lwd=2)
segments(min(data3\$qc),max(data3\$number)*1.1,min(data3\$qc),40,lwd=2)
axis(2,col="grey40",col.axis="grey40",las=2,pos=min(data3\$qc),lwd=2)
axis(1,col="grey40",col.axis="grey40",pos=0,lwd=2,at=xlabel,labels=xlabel)
mtext(side=2,"Number of Reads",line=3,cex=1.5)
mtext(side=3,"Per sequence quality scores",line=1,cex=1.8)
mtext(side=1,"Mean Sequence Quality(Phred Score)",line=2,cex=1.5)
################3.Per_sequence_quality_scores: end

par(mfrow = c(1,1))
data4 <- read.table("$outdir/result/read1.Per_base_sequence_content.xls",head=T,sep = "\\t")
xlabel <- c(1,seq(10,length(row.names(data4)),10))
ylabel <- c(seq(0,max(data4[,-1]),10))
par(mar=c(4,4,4,5),xpd=TRUE)
plot(data4\$pos,data4\$A,ylab="",xlab="",cex.lab=1.5,yaxt="n",xaxt="n",axes=F,main="",ylim=c(0,max(data4[,-1])),xlim=c(0,max(data4\$pos)),cex=0)
axis(2,col="grey70",col.axis="grey40",las=2,pos=1,lwd=2,at = ylabel,labels = ylabel)
axis(1,col="grey70",col.axis="grey40",pos=0,lwd=2,at = xlabel,labels = xlabel)
segments(1,0,1,max(data4[,-1]),lwd=2,col="grey70")
mtext(side=2,"Percent base content",line=2,cex=1.5)
mtext(side=1,"Position in read(bp)",line=2,cex=1.5)
mtext(side=3,"Per base sequence content",line=1,cex=1.8)
text(length(row.names(data4))*1.07,(max(data4[,-1]))/2+(max(data4[,-1]))*0.2,labels="%A",bty = "n",col="#dc143c",cex=1.2)
text(length(row.names(data4))*1.07,(max(data4[,-1]))/2+(max(data4[,-1]))*0.1,labels="%T",bty = "n",col="#1e90ff",cex=1.2)
text(length(row.names(data4))*1.07,(max(data4[,-1]))/2+(max(data4[,-1]))*0.0,labels="%C",bty = "n",col="#f4a460",cex=1.2)
text(length(row.names(data4))*1.07,(max(data4[,-1]))/2-(max(data4[,-1]))*0.1,labels="%G",bty = "n",col="#32cd32",cex=1.2)
text(length(row.names(data4))*1.07,(max(data4[,-1]))/2-(max(data4[,-1]))*0.2,labels="%N",bty = "n",col="#696969",cex=1.2)
par(new=T)
par(mar=c(4,4,4,5),xpd=TRUE)
plot(data4\$pos,data4\$A,ylab="",xlab="",cex.lab=1.5,yaxt="n",xaxt="n",axes=F,main="",ylim=c(0,max(data4[,-1])),xlim=c(0,max(data4\$pos)),cex=0)
lines(data4\$pos,data4\$A,lty=1,col="#dc143c",lwd=2)
lines(data4\$pos,data4\$T,lty=1,col="#1e90ff",lwd=2)
lines(data4\$pos,data4\$C,lty=1,col="#f4a460",lwd=2)
lines(data4\$pos,data4\$G,lty=1,col="#32cd32",lwd=2)
lines(data4\$pos,data4\$N,lty=1,col="#696969",lwd=2)
#########4.Per base sequence content :end     

par(mfrow = c(1,1))
data5 <- read.table("$outdir/result/read1.Per_sequence_GC_content.xls",head=T,sep = "\\t")
xlabel <- c(0,seq(10,100,10))
par(mar=c(4,5,4,3),xpd=TRUE)
plot(data5\$gc,data5\$number,ylab="",xlab="",xaxt="n",ylim=c(0,max(data5\$number)*1.2),axes=F,cex=0)
lines(data5\$gc,data5\$number,lty=1,col="#dc143c",lwd=2)
par(new=T)
par(mar=c(4,5,4,3),xpd=TRUE)
axis(2,col="grey40",col.axis="grey40",las=2,pos=min(data5\$gc),lwd=2)
axis(1,col="grey40",col.axis="grey40",pos=0.01,lwd=2,at=xlabel,label=xlabel)
segments(0,0,0,max(data5\$number)*1.1,lwd=2)
mtext(side=2,"Number of Reads",line=3,cex=1.5)
mtext(side=3,"Per sequence GC content",line=0,cex=1.8)
mtext(side=1,"GC Content(%)",line=2,cex=1.5)
########5.Per_sequence_GC_content :end
dev.off()
CMD

system("$Rscript $outdir/result/SE.QC_Statistics.R 2>/dev/null");
system("$Convert -density 300 -resize 30% $outdir/result/SE.QC_Statistics.pdf $outdir/result/SE.QC_Statistics.png");
}

else
{
open my $fh_rcode,">","$outdir/result/PE.QC_Statistics.R" or die $!;
print $fh_rcode <<CMD;
pdf("$outdir/result/PE.QC_Statistics.pdf",8,5)
par(mfrow = c(1,1))
data1 <- read.table("$outdir/result/read1.Basic_Statistics.xls",head=T,sep="\\t")
len <- length(row.names(data1))
xc <- seq(0,len+2,1)
yc <- seq(0,len+2,1)
plot(xc,yc,ylab="",xlab="",cex.lab=1.5,yaxt="n",xaxt="n",axes=F,main="",cex=0)
cc1 <- seq(0.5,len+0.5,0.05)
cc2 <- seq(len+0.5,len+1.5,0.05)
for (i in cc1)
{
        segments(0,i,len+2,i,lwd=4,col="#DEEBF7")
}
for (i in cc2)
{
        segments(0,i,len+2,i,lwd=4,col="#9ECAE1")
}
segments(0,len+0.5,len+2,len+0.5,lwd=3,col="grey20")
segments(0,len+1.5,len+2,len+1.5,lwd=3,col="grey20")
segments(0,0.5,len+2,0.5,lwd=3,col="grey20")
mtext(side=3,"Basic Statistics(read1)",line=0,cex=1.8)
text(3,length(row.names(data1))+1,labels="Measure")
text(7,length(row.names(data1))+1,labels="Value")
for (i in 1:length(row.names(data1)))
{
        text(3,(length(row.names(data1))+1- i),labels=data1\$Measure[i])
        text(7,(length(row.names(data1))+1- i),labels=data1\$Value[i])
}

par(mfrow = c(1,1))
data1 <- read.table("$outdir/result/read2.Basic_Statistics.xls",head=T,sep="\\t")
len <- length(row.names(data1))
xc <- seq(0,len+2,1)
yc <- seq(0,len+2,1)
plot(xc,yc,ylab="",xlab="",cex.lab=1.5,yaxt="n",xaxt="n",axes=F,main="",cex=0)
cc1 <- seq(0.5,len+0.5,0.05)
cc2 <- seq(len+0.5,len+1.5,0.05)
for (i in cc1)
{
        segments(0,i,len+2,i,lwd=4,col="#DEEBF7")
}
for (i in cc2)
{
        segments(0,i,len+2,i,lwd=4,col="#9ECAE1")
}
segments(0,len+0.5,len+2,len+0.5,lwd=3,col="grey20")
segments(0,len+1.5,len+2,len+1.5,lwd=3,col="grey20")
segments(0,0.5,len+2,0.5,lwd=3,col="grey20")
mtext(side=3,"Basic Statistics(read2)",line=0,cex=1.8)
text(3,length(row.names(data1))+1,labels="Measure")
text(7,length(row.names(data1))+1,labels="Value")
for (i in 1:length(row.names(data1)))
{
        text(3,(length(row.names(data1))+1- i),labels=data1\$Measure[i])
        text(7,(length(row.names(data1))+1- i),labels=data1\$Value[i])
}

#########1.Basic_Statistics :end

par(mfrow = c(1,1))
data2 <- read.table("$outdir/result/read1.Per_base_sequence_quality.xls",head=T,sep="\\t")
pos_len <- length(row.names(data2))
k <- seq(0,40-40/pos_len,40/pos_len)
w0 <- seq(0,20,0.05)
w20 <- seq(20,30,0.05)
w30 <- seq(30,40,0.05)
par(mar=c(4,5,4,2),xpd=TRUE)
plot(data2\$Position,k,xlab="",ylab="",xaxt="n",yaxt="n",cex=0,ylim = c(0,40),axes = FALSE)
for (i in 1:length(w0))
{
        segments(0,w0[i],max(data2\$Position)+1,w0[i],lwd=2,col=adjustcolor("#C6DBEF",alpha.f = 0.2))
}
for (i in 1:length(w20))
{
        segments(0,w20[i],max(data2\$Position)+1,w20[i],lwd=2,col=adjustcolor("#6BAED6",alpha.f = 0.2))
}
for (i in 1:length(w30))
{
        segments(0,w30[i],max(data2\$Position)+1,w30[i],lwd=2,col=adjustcolor("#2171B5",alpha.f = 0.2))
}
par(new=T)
par(mar=c(4,5,4,2),xpd=TRUE)
xlabel <-c(1,seq(10,max(data2\$Position),10))
ylabel <- seq(0,40,10)
seg.back <- seq(-0.3,0.3,0.02)
plot(data2\$Position,k,xlab="",ylab="",xaxt="n",yaxt="n",cex=0,ylim = c(0,40),axes = FALSE)
segments(0,0,0,40,lwd=2)
segments(0,0,max(data2\$Position)+1,0,lwd=2)
segments(0,20,max(data2\$Position)+1,20,lwd=2,col="grey70",lty=2)
segments(0,30,max(data2\$Position)+1,30,lwd=2,col="grey70",lty=2)
for (i in 1:pos_len)
{
        segments(i-0.2,data2\$X10th.Percentile[i],i+0.2,data2\$X10th.Percentile[i],lwd=1,col="grey50")
        segments(i-0.2,data2\$X90th.Percentile[i],i+0.2,data2\$X90th.Percentile[i],lwd=1,col="grey50")
        segments(i,data2\$X10th.Percentile[i],i,data2\$X90th.Percentile[i],lwd=1,col="grey50",lty=2)
	for (j in seg.back)
        {
                segments(i+j,data2\$Lower.Quartile[i],i+j,data2\$Upper.Quartile[i],lwd=0.5,col="#FB9A99")
        }
        segments(i-0.3,data2\$Median[i],i+0.3,data2\$Median[i],lwd=1,col=adjustcolor("#E41A1C",alpha.f = 0.9))
}  
lines(data2\$Mean, lty =1, col = "#6A51A3")
axis(1,lwd=2,col="grey40",col.axis="grey40",pos=0,at=xlabel,labels =xlabel)
axis(2,lwd=2,col="grey40",col.axis="grey40",pos=0,at=ylabel,label=ylabel,las=2)
mtext(side=2,"Quality Scores",line=3,cex=1.5)
mtext(side=3,"Per base sequence quality(read1)",line=1,cex=1.8)
mtext(side=1,"Position in read(bp)",line=2,cex=1.5)

par(mfrow = c(1,1))
data2 <- read.table("$outdir/result/read2.Per_base_sequence_quality.xls",head=T,sep="\\t")
pos_len <- length(row.names(data2))
k <- seq(0,40-40/pos_len,40/pos_len)
w0 <- seq(0,20,0.05)
w20 <- seq(20,30,0.05)
w30 <- seq(30,40,0.05)
par(mar=c(4,5,4,2),xpd=TRUE)
plot(data2\$Position,k,xlab="",ylab="",xaxt="n",yaxt="n",cex=0,ylim = c(0,40),axes = FALSE)
for (i in 1:length(w0))
{
        segments(0,w0[i],max(data2\$Position)+1,w0[i],lwd=2,col=adjustcolor("#C6DBEF",alpha.f = 0.2))
}
for (i in 1:length(w20))
{
        segments(0,w20[i],max(data2\$Position)+1,w20[i],lwd=2,col=adjustcolor("#6BAED6",alpha.f = 0.2))
}
for (i in 1:length(w30))
{
        segments(0,w30[i],max(data2\$Position)+1,w30[i],lwd=2,col=adjustcolor("#2171B5",alpha.f = 0.2))
}
par(new=T)
par(mar=c(4,5,4,2),xpd=TRUE)
xlabel <-c(1,seq(10,max(data2\$Position),10))
ylabel <- seq(0,40,10)
seg.back <- seq(-0.3,0.3,0.02)
plot(data2\$Position,k,xlab="",ylab="",xaxt="n",yaxt="n",cex=0,ylim = c(0,40),axes = FALSE)
segments(0,0,0,40,lwd=2)
segments(0,0,max(data2\$Position)+1,0,lwd=2)
segments(0,20,max(data2\$Position)+1,20,lwd=2,col="grey70",lty=2)
segments(0,30,max(data2\$Position)+1,30,lwd=2,col="grey70",lty=2)
for (i in 1:pos_len)
{
        segments(i-0.2,data2\$X10th.Percentile[i],i+0.2,data2\$X10th.Percentile[i],lwd=1,col="grey50")
        segments(i-0.2,data2\$X90th.Percentile[i],i+0.2,data2\$X90th.Percentile[i],lwd=1,col="grey50")
        segments(i,data2\$X10th.Percentile[i],i,data2\$X90th.Percentile[i],lwd=1,col="grey50",lty=2)
        for (j in seg.back)
        {
                segments(i+j,data2\$Lower.Quartile[i],i+j,data2\$Upper.Quartile[i],lwd=0.5,col="#FB9A99")
        }
        segments(i-0.3,data2\$Median[i],i+0.3,data2\$Median[i],lwd=1,col=adjustcolor("#E41A1C",alpha.f = 0.9))
}  
lines(data2\$Mean, lty =1, col = "#6A51A3")
axis(1,lwd=2,col="grey40",col.axis="grey40",pos=0,at=xlabel,labels =xlabel)
axis(2,lwd=2,col="grey40",col.axis="grey40",pos=0,at=ylabel,label=ylabel,las=2)
mtext(side=2,"Quality Scores",line=3,cex=1.5)
mtext(side=3,"Per base sequence quality(read2)",line=1,cex=1.8)
mtext(side=1,"Position in read(bp)",line=2,cex=1.5)

#################2.Per_sequence_quality_scores :end

par(mfrow = c(1,1))
data3 <- read.table("$outdir/result/read1.Per_sequence_quality_scores.xls",head=T,sep = "\\t")
xlabel <- c(min(data3\$qc),seq(10,40,5))
par(mar=c(4,5,6,4),xpd=TRUE)
plot(data3\$qc,data3\$number,ylab="",xlab="",xaxt="n",ylim=c(0,max(data3\$number)*1.2),xlim=c(min(data3\$qc),38),axes=F,cex=0)
lines(data3\$qc,data3\$number,lty=1,col="#dc143c",lwd=2)
par(new=T)
par(mar=c(4,5,6,4),xpd=TRUE)
segments(min(data3\$qc),0,max(data3\$qc)+1,0,lwd=2)
segments(min(data3\$qc),max(data3\$number)*1.1,min(data3\$qc),40,lwd=2)
axis(2,col="grey40",col.axis="grey40",las=2,pos=min(data3\$qc),lwd=2)
axis(1,col="grey40",col.axis="grey40",pos=0,lwd=2,at=xlabel,labels=xlabel)
mtext(side=2,"Number of Reads",line=3,cex=1.5)
mtext(side=3,"Per sequence quality scores(read1)",line=1,cex=1.8)
mtext(side=1,"Mean Sequence Quality(Phred Score)",line=2,cex=1.5)

par(mfrow = c(1,1))
data3 <- read.table("$outdir/result/read2.Per_sequence_quality_scores.xls",head=T,sep = "\\t")
xlabel <- c(min(data3\$qc),seq(10,40,5))
par(mar=c(4,5,6,4),xpd=TRUE)
plot(data3\$qc,data3\$number,ylab="",xlab="",xaxt="n",ylim=c(0,max(data3\$number)*1.2),xlim=c(min(data3\$qc),38),axes=F,cex=0)
lines(data3\$qc,data3\$number,lty=1,col="#dc143c",lwd=2)
par(new=T)
par(mar=c(4,5,6,4),xpd=TRUE)
segments(min(data3\$qc),0,max(data3\$qc)+1,0,lwd=2)
segments(min(data3\$qc),max(data3\$number)*1.1,min(data3\$qc),40,lwd=2)
axis(2,col="grey40",col.axis="grey40",las=2,pos=min(data3\$qc),lwd=2)
axis(1,col="grey40",col.axis="grey40",pos=0,lwd=2,at=xlabel,labels=xlabel)
mtext(side=2,"Number of Reads",line=3,cex=1.5)
mtext(side=3,"Per sequence quality scores(read2)",line=1,cex=1.8)
mtext(side=1,"Mean Sequence Quality(Phred Score)",line=2,cex=1.5)

################3.Per_sequence_quality_scores: end

par(mfrow = c(1,1))
data4 <- read.table("$outdir/result/read1.Per_base_sequence_content.xls",head=T,sep = "\\t")
xlabel <- c(1,seq(10,length(row.names(data4)),10))
ylabel <- c(seq(0,max(data4[,-1]),10))
par(mar=c(4,4,4,5),xpd=TRUE)
plot(data4\$pos,data4\$A,ylab="",xlab="",cex.lab=1.5,yaxt="n",xaxt="n",axes=F,main="",ylim=c(0,max(data4[,-1])),xlim=c(0,max(data4\$pos)),cex=0)
axis(2,col="grey70",col.axis="grey40",las=2,pos=1,lwd=2,at = ylabel,labels = ylabel)
axis(1,col="grey70",col.axis="grey40",pos=0,lwd=2,at = xlabel,labels = xlabel)
segments(1,0,1,max(data4[,-1]),lwd=2,col="grey70")
mtext(side=2,"Percent base content",line=2,cex=1.5)
mtext(side=1,"Position in read(bp)",line=2,cex=1.5)
mtext(side=3,"Per base sequence content(read1)",line=1,cex=1.8)
text(length(row.names(data4))*1.07,(max(data4[,-1]))/2+(max(data4[,-1]))*0.2,labels="%A",bty = "n",col="#dc143c",cex=1.2)
text(length(row.names(data4))*1.07,(max(data4[,-1]))/2+(max(data4[,-1]))*0.1,labels="%T",bty = "n",col="#1e90ff",cex=1.2)
text(length(row.names(data4))*1.07,(max(data4[,-1]))/2+(max(data4[,-1]))*0.0,labels="%C",bty = "n",col="#f4a460",cex=1.2)
text(length(row.names(data4))*1.07,(max(data4[,-1]))/2-(max(data4[,-1]))*0.1,labels="%G",bty = "n",col="#32cd32",cex=1.2)
text(length(row.names(data4))*1.07,(max(data4[,-1]))/2-(max(data4[,-1]))*0.2,labels="%N",bty = "n",col="#696969",cex=1.2)
par(new=T)
par(mar=c(4,4,4,5),xpd=TRUE)
plot(data4\$pos,data4\$A,ylab="",xlab="",cex.lab=1.5,yaxt="n",xaxt="n",axes=F,main="",ylim=c(0,max(data4[,-1])),xlim=c(0,max(data4\$pos)),cex=0)
lines(data4\$pos,data4\$A,lty=1,col="#dc143c",lwd=2)
lines(data4\$pos,data4\$T,lty=1,col="#1e90ff",lwd=2)
lines(data4\$pos,data4\$C,lty=1,col="#f4a460",lwd=2)
lines(data4\$pos,data4\$G,lty=1,col="#32cd32",lwd=2)
lines(data4\$pos,data4\$N,lty=1,col="#696969",lwd=2)

par(mfrow = c(1,1))
data4 <- read.table("$outdir/result/read2.Per_base_sequence_content.xls",head=T,sep = "\\t")
xlabel <- c(1,seq(10,length(row.names(data4)),10))
ylabel <- c(seq(0,max(data4[,-1]),10))
par(mar=c(4,4,4,5),xpd=TRUE)
plot(data4\$pos,data4\$A,ylab="",xlab="",cex.lab=1.5,yaxt="n",xaxt="n",axes=F,main="",ylim=c(0,max(data4[,-1])),xlim=c(0,max(data4\$pos)),cex=0)
axis(2,col="grey70",col.axis="grey40",las=2,pos=1,lwd=2,at = ylabel,labels = ylabel)
axis(1,col="grey70",col.axis="grey40",pos=0,lwd=2,at = xlabel,labels = xlabel)
segments(1,0,1,max(data4[,-1]),lwd=2,col="grey70")
mtext(side=2,"Percent base content",line=2,cex=1.5)
mtext(side=1,"Position in read(bp)",line=2,cex=1.5)
mtext(side=3,"Per base sequence content(read2)",line=1,cex=1.8)
text(length(row.names(data4))*1.07,(max(data4[,-1]))/2+(max(data4[,-1]))*0.2,labels="%A",bty = "n",col="#dc143c",cex=1.2)
text(length(row.names(data4))*1.07,(max(data4[,-1]))/2+(max(data4[,-1]))*0.1,labels="%T",bty = "n",col="#1e90ff",cex=1.2)
text(length(row.names(data4))*1.07,(max(data4[,-1]))/2+(max(data4[,-1]))*0.0,labels="%C",bty = "n",col="#f4a460",cex=1.2)
text(length(row.names(data4))*1.07,(max(data4[,-1]))/2-(max(data4[,-1]))*0.1,labels="%G",bty = "n",col="#32cd32",cex=1.2)
text(length(row.names(data4))*1.07,(max(data4[,-1]))/2-(max(data4[,-1]))*0.2,labels="%N",bty = "n",col="#696969",cex=1.2)
par(new=T)
par(mar=c(4,4,4,5),xpd=TRUE)
plot(data4\$pos,data4\$A,ylab="",xlab="",cex.lab=1.5,yaxt="n",xaxt="n",axes=F,main="",ylim=c(0,max(data4[,-1])),xlim=c(0,max(data4\$pos)),cex=0)
lines(data4\$pos,data4\$A,lty=1,col="#dc143c",lwd=2)
lines(data4\$pos,data4\$T,lty=1,col="#1e90ff",lwd=2)
lines(data4\$pos,data4\$C,lty=1,col="#f4a460",lwd=2)
lines(data4\$pos,data4\$G,lty=1,col="#32cd32",lwd=2)
lines(data4\$pos,data4\$N,lty=1,col="#696969",lwd=2)

#########4.Per base sequence content :end

par(mfrow = c(1,1))
data5 <- read.table("$outdir/result/read1.Per_sequence_GC_content.xls",head=T,sep = "\\t")
xlabel <- c(0,seq(10,100,10))
par(mar=c(4,5,5,4),xpd=TRUE)
plot(data5\$gc,data5\$number,ylab="",xlab="",xaxt="n",ylim=c(0,max(data5\$number)*1.2),axes=F,cex=0)
lines(data5\$gc,data5\$number,lty=1,col="#dc143c",lwd=2)
data6 <- read.table("$outdir/result/read2.Per_sequence_GC_content.xls",head=T,sep = "\\t")
lines(data6\$gc,data6\$number,lty=1,col="#1e90ff",lwd=2)
par(new=T)
par(mar=c(4,5,5,4),xpd=TRUE)
axis(2,col="grey40",col.axis="grey40",las=2,pos=min(data5\$gc),lwd=2)
axis(1,col="grey40",col.axis="grey40",pos=0.01,lwd=2,at=xlabel,label=xlabel)
segments(0,0,0,max(data5\$number)*1.1,lwd=2)
mtext(side=2,"Number of Reads",line=3,cex=1.5)
mtext(side=3,"Per sequence GC content",line=0,cex=1.8)
mtext(side=1,"GC Content(%)",line=2,cex=1.5)
text(max(data5\$gc)*1.07,(max(data5\$number))/2+(max(data5\$number))*0.3,labels="read1",bty = "n",col="#dc143c",cex=1.4)
text(max(data5\$gc)*1.07,(max(data5\$number))/2+(max(data5\$number))*0.1,labels="read2",bty = "n",col="#1e90ff",cex=1.4)
########5.Per_sequence_GC_content :end
dev.off()
CMD

system("$Rscript $outdir/result/PE.QC_Statistics.R 2>/dev/null");
system("$Convert -density 300 -resize 30% $outdir/result/PE.QC_Statistics.pdf $outdir/result/PE.QC_Statistics.png");
}
system("echo end plot: `date`");



sub start_thread
{
	my ($N_per,$GC_per,$Q20,$Q30,$N_num,$GC,%g_num,%a_num,%c_num,%t_num,%n_num);
	$N_num = 0;
	$GC = 0;
	my $file_seg = shift;
	my $read = shift;
	my $seq_number = shift;
	open (FQ1,"$file_seg")or die $!;
	while(<FQ1>)
	{
		$seq_gc = 0;
		$base_line=<FQ1>;
		<FQ1>;
		$qual_line=<FQ1>;
		my @tmp = split /\n/;
		my @qual = split //,$qual_line;pop(@qual);
		my @base = split //,$base_line;pop(@base);
		foreach $pos (1..$read_length)
		{
			if ($base[$pos-1]=~m/G/){$g_num{$pos} ++;$GC ++;$seq_gc ++;}
			if ($base[$pos-1]=~m/C/){$c_num{$pos} ++;$GC ++;$seq_gc ++;}
			if ($base[$pos-1]=~m/A/){$a_num{$pos} ++;}
			if ($base[$pos-1]=~m/T/){$t_num{$pos} ++;}
			if ($base[$pos-1]=~m/N/){$n_num{$pos} ++;$N_num ++;}
			my $q = ord $qual[$pos-1];
                        $q = $q - $asi_qc;
                        $qual_count_pos{$pos}{$q} ++;
                        $qual_count{$q} ++;
		}
		$sequence_gc{$seq_gc} ++;
	}
    	$Q20 = 0;
    	$Q30 = 0;
    	foreach $qc (20..40)
   	{
        	$Q20 = $Q20 + $qual_count{$qc};
    	}
    	foreach $qc (30..40)
    	{
    		$Q30 = $Q30 + $qual_count{$qc};
    	}
	$Q20 = sprintf ("%.2f",($Q20/($seq_number*$read_length))*100);
   	$Q30 = sprintf ("%.2f",($Q30/($seq_number*$read_length))*100);
	$N_per = sprintf ("%.2f",($N_num/($seq_number*$read_length))*100);	
	$GC_per = sprintf ("%.2f",($GC/($seq_number*$read_length))*100);
    	if ($read eq "read1") {print POUT1 "$GC_per\t$N_per\t$Q20\t$Q30\n";}
	if ($read eq "read2") {print ROUT1 "$GC_per\t$N_per\t$Q20\t$Q30\n";}
	my(%pos_qual_ranking,%pos_qual_mean,%pos_qual_median,%pos_qual_lower_quartile,%pos_qual_upper_quartile,%pos_qual_10th_percentile,%pos_qual_90th_percentile);
	foreach $pos (1..$read_length)
	{
		my %pos_qual_rangking = ();
		my $qual_ranking = 0;
		my $qual_total = 0;
		foreach $qual (0..40)
		{
			if(exists $qual_count_pos{$pos}{$qual}){$q_num = $qual_count_pos{$pos}{$qual};}
			else {$q_num=0;}
			$qual_total = $qual_total + $q_num * $qual;
			$qual_ranking = $qual_ranking + $q_num;
			$pos_qual_ranking{$qual} = $qual_ranking;
		}
		my $qual_median = $qual_ranking *0.5;
		my $qual_lower_quartile = $qual_ranking *0.25;
		my $qual_upper_quartile = $qual_ranking *0.75;
		my $qual_10th_percentile = $qual_ranking *0.1;
		my $qual_90th_percentile = $qual_ranking *0.9;
		my $mean_qual = $qual_total/$qual_ranking;
		$pos_qual_mean{$pos} = $mean_qual;
		foreach $qual(0..40)
		{
			$qual_ranking = $pos_qual_ranking{$qual};
			$lower_qual = $qual -1;
			$smaller_qual_ranking = $pos_qual_ranking{$lower_qual};
			if ($qual_ranking >= $qual_median and $smaller_qual_ranking < $qual_median)
			{
				$pos_qual_median{$pos} = $qual;
			}
			if ($qual_ranking >= $qual_lower_quartile and $smaller_qual_ranking < $qual_lower_quartile)
			{
				$pos_qual_lower_quartile{$pos} = $qual;
			}
			if ($qual_ranking >= $qual_upper_quartile and $smaller_qual_ranking < $qual_upper_quartile)
			{
				$pos_qual_upper_quartile{$pos} = $qual;
			}
			if ($qual_ranking >= $qual_10th_percentile and $smaller_qual_ranking < $qual_10th_percentile)
			{
				$pos_qual_10th_percentile{$pos} = $qual;
			}
			if ($qual_ranking >= $qual_90th_percentile and $smaller_qual_ranking < $qual_90th_percentile)
			{
				$pos_qual_90th_percentile{$pos} = $qual;
			}
		}
	}
	if ($read eq "read1")
	{
	foreach $pos (1..$read_length)
	{
        	print POUT2 "$pos\tMean\t$pos_qual_mean{$pos}\n";
		print POUT2 "$pos\tMedian\t$pos_qual_median{$pos}\n";
		print POUT2 "$pos\tLower Quartile\t$pos_qual_lower_quartile{$pos}\n";
		print POUT2 "$pos\tUpper Quartile\t$pos_qual_upper_quartile{$pos}\n";
		print POUT2 "$pos\t10th Percentile\t$pos_qual_10th_percentile{$pos}\n";
		print POUT2 "$pos\t90th Percentile\t$pos_qual_90th_percentile{$pos}\n";
	}
	foreach $q (sort {$a <=> $b} keys %qual_count)
	{
        	print POUT3 "$q\t$qual_count{$q}\n";
	}
	foreach $pos (sort {$a <=> $b} keys %g_num)
	{
		if (!exists $n_num{$pos}){$n_num{$pos} =0;}
		$a = ($a_num{$pos}/$seq_number)*100;
		$t = ($t_num{$pos}/$seq_number)*100;
		$c = ($c_num{$pos}/$seq_number)*100;
		$g = ($g_num{$pos}/$seq_number)*100;
		$n = ($n_num{$pos}/$seq_number)*100;
		print POUT4 "$pos\tA\t$a\n";
		print POUT4 "$pos\tT\t$t\n";
		print POUT4 "$pos\tC\t$c\n";
		print POUT4 "$pos\tG\t$g\n";
		print POUT4 "$pos\tN\t$n\n";
	}
	foreach $seq_gc (sort {$a <=> $b} keys %sequence_gc)
	{
        	$content = ($seq_gc/$read_length)*100;
        	print POUT5 "$content\t$sequence_gc{$seq_gc}\n";
	}
	}
	
	if ($read eq "read2")
        {
        foreach $pos (1..$read_length)
        {
                print ROUT2 "$pos\tMean\t$pos_qual_mean{$pos}\n";
                print ROUT2 "$pos\tMedian\t$pos_qual_median{$pos}\n";
                print ROUT2 "$pos\tLower Quartile\t$pos_qual_lower_quartile{$pos}\n";
                print ROUT2 "$pos\tUpper Quartile\t$pos_qual_upper_quartile{$pos}\n";
                print ROUT2 "$pos\t10th Percentile\t$pos_qual_10th_percentile{$pos}\n";
                print ROUT2 "$pos\t90th Percentile\t$pos_qual_90th_percentile{$pos}\n";
        }
        foreach $q (sort {$a <=> $b} keys %qual_count)
        {
                print ROUT3 "$q\t$qual_count{$q}\n";
        }
        foreach $pos (sort {$a <=> $b} keys %g_num)
        {
                if (!exists $n_num{$pos}){$n_num{$pos} =0;}
                $a = ($a_num{$pos}/$seq_number)*100;
                $t = ($t_num{$pos}/$seq_number)*100;
                $c = ($c_num{$pos}/$seq_number)*100;
                $g = ($g_num{$pos}/$seq_number)*100;
                $n = ($n_num{$pos}/$seq_number)*100;
                print ROUT4 "$pos\tA\t$a\n";
                print ROUT4 "$pos\tT\t$t\n";
                print ROUT4 "$pos\tC\t$c\n";
                print ROUT4 "$pos\tG\t$g\n";
                print ROUT4 "$pos\tN\t$n\n";
        }
        foreach $seq_gc (sort {$a <=> $b} keys %sequence_gc)
        {
                $content = ($seq_gc/$read_length)*100;
                print ROUT5 "$content\t$sequence_gc{$seq_gc}\n";
        }
        }
}
if(!defined $file2){system("rm $outdir/result/read1.fq");}
else {system("rm $outdir/result/read1.fq & rm $outdir/result/read2.fq");}
system("rm -rf $outdir/process");
