package panzhong_win7;
use Exporter;
use File::Copy;
use Time::Local;
use Cwd;
my $mydir = getcwd();
our @ISA = qw(Exporter);
our @EXPORT = qw(get_sub_list creat_database get_table_from_sql get_table_with_head drop_table);
push @EXPORT, qw(input_table_auto_head input_table_sql input_poss_inf_sql input_poss_sql check_sql_log);
push @EXPORT, qw(input_bed input_psl input_blast input_novoalign_native input_sam input_gbCdnaInfo sqlserver_table_exists sqlserver_table_line_count);
push @EXPORT, qw(fasta_out fa_out anti_reverse circRNA_forTargetScan_out get_targetscan_format get_sub_fasta get_circRNA_sequence_forTargetScan get_fasta_length);
push @EXPORT, qw(gettime get_day time_interval);
push @EXPORT, qw(two_excel_win32 add_column read_config txt2xlsx_win32 get_taxid get_expression_all_txt get_md5_files get_differentially_expressed_circRNA_list_win32);
push @EXPORT, qw(txt2xlsx get_differentially_expressed_circRNA_list two_excel_writer);


use Digest::MD5;
use Spreadsheet::ParseExcel;
use Excel::Writer::XLSX;


use Win32::OLE;
#use Win32::OLE::Const 'Microsoft Excel';
my $server='YXSY11\SQL2012EXPRESS';
my $user ='sa';
my $pw='110883';
my $shortdir="d:";

##################       output
=pod
              &get_table_from_sql($database,$table);
              &get_table_with_head($database,$table);
=cut
##################       input
=pod
            my @columns=qw(chrom txStart txEnd);
            my @types= split(" ","varchar(255) int(10) int(10)");
               &input_table_sql($database,$filename,$table,\@columns,\@types);
               &input_table_by_poss($database,$file,$table);
               &input_table_auto_head($database,$file,$table);
                          #######
                          &input_poss_inf_sql($database,$file,$table);
                          &input_poss_sql($database,$file,$table);
                          &input_bed($database,$filename,$table,6);
                          &input_psl($database,$psl,$table);
                          &input_blast($database,$filename,$table);
                          &input_novoalign_native($database,$file,$table);
                          &input_sam($database,$file,$table);
                          &input_gbCdnaInfo($database,$file,$table);
=cut

sub get_circRNA_sequence_forTargetScan
{
open (LIST, "$_[0].txt") or die "error(input1):$!";
open(TEMP, ">$_[0]_forTargetScan.fa") or die "error (output1):$!";
open(MISS, ">$_[0].mis") or die "error (output2):$!";

my %sublist = ();
my $line;
while ( $line= <LIST>)
{     chomp $line;
      my @name = split(/\t/,$line);
      my $name = $name[0];
     # print $name,"\n";
      $sublist{$name} = 1;
}
my $count = keys(%sublist);
# print "there are $count unique record in the list file!\n";

my %gotlist = ();
my $name;
foreach my $source(@{$_[1]})
{
open (FASTA, "$source") or die "error(input2):$!";
my $seq;
OUTER: while ( my $fastaline= <FASTA> )
{
	$name="KKKKKKK";
		if($fastaline =~ />/)
	{         $fastaline =~ />([\w\:\|\_\.\-\+\d]+)/;
                  $name = $1;
                  $seq = '';
                 if ( (exists($sublist{$name}))  and (not exists($gotlist{$name})))
                 {

	            	  while($fastaline= <FASTA>)
	           	    {
                                    if(not $fastaline =~ />/)
	               		   {
                                    $seq .= $fastaline;
                                    }
	          		     else
                                    {
                                     $seq = &circRNA_forTargetScan_out($seq);
                                    print TEMP ">$name\n$seq";
                                    $gotlist{$name} = 1;
                                    $seq = '';
                                    redo OUTER;
                                    }
         		    }#while
	         }#if
        }##if
}#while OUT
                          if($name ne "KKKKKKK")
                          {
                            if ( (exists($sublist{$name}))  and (not exists($gotlist{$name})))
                               {
                                  $seq = &circRNA_forTargetScan_out($seq);
                                  print TEMP ">$name\n$seq";
                                  $gotlist{$name} = 1;
                                 $seq = '';
                               }
                           }
close FASTA
}#foreach

$count = keys %gotlist;
          print "there are $count CircRNA sequences transformed into TargetScan format!\n";
foreach(keys %sublist)
{
     my $tag = $_;
       if(not exists($gotlist{$tag}))
         {
           print MISS "$tag","\n";
         }#if
}#foreach
 close (LIST);
 close (TEMP);
 close (MISS);
}

sub get_sub_fasta
{
open (LIST, "$_[0].txt") or die "error(input1):$!";
open(TEMP, ">$_[0].fa") or die "error (output1):$!";
open(MISS, ">$_[0].mis") or die "error (output2):$!";

my %sublist = ();
my $line;
while ( $line= <LIST>)
{     chomp $line;
      my @name = split(/\t/,$line);
      my $name = $name[$_[2]];
     # print $name,"\n";
      $sublist{$name} = 1;
}
my $count = keys(%sublist);
print "there are $count unique record in the list file!\n";

my %gotlist = ();
my @input;
if(ref($_[1]))
{
	@input = @{$_[1]};
	}
else
{
	push @input,$_[1];
	}

foreach my $source(@input)
{
open (FASTA, "$source") or die "error(input2):$!";
my $name;
my $seq;
OUTER: while ( my $fastaline= <FASTA> )
{
	$name="KKKKKKK";
		if($fastaline =~ />/)
	{         $fastaline =~ />([\w\:\_\.\-\+\d]+)/;
                  $name = $1;
                  $seq = '';
                 if ( (exists($sublist{$name}))  and (not exists($gotlist{$name})))
                 {

	            	  while($fastaline= <FASTA>)
	           	    {
                                    if(not $fastaline =~ />/)
	               		   {
                                    $seq .= $fastaline;
                                    }
	          		   else
                                    {
                                    $seq = &fasta_out($seq);
                                    print TEMP ">$name\n$seq";
                                    $gotlist{$name} = 1;
                                    $seq = '';
                                    redo OUTER;
                                    }
         		    }#while
	         }#if
        }##if
}#while OUT
                          if($name ne "KKKKKKK")
                          {
                            if ( (exists($sublist{$name}))  and (not exists($gotlist{$name})))
                               {
                                  $seq = &fasta_out($seq);
                                  print TEMP ">$name\n$seq";
                                  $gotlist{$name} = 1;
                                 $seq = '';
                               }
                           }
close FASTA
}#foreach

$count = keys %gotlist;
print "there are $count unique record in the gotlist file!\n";
foreach(keys %sublist)
{
     my $tag = $_;
       if(not exists($gotlist{$tag}))
         {
           print MISS "$tag","\n";
         }#if
}#foreach
 close (LIST);
 close (TEMP);
 close (MISS);
}


sub get_fasta_length
{
print "inference the $_[0].fa sequence length..................!\n";
open (SEQUENCE1, "$_[0].fa") or die "error(input):$!";
open(SEQUENCE2, ">temp.fa") or die "error (output):$!";
while(my $seq = <SEQUENCE1>)
{    if( $seq =~ />/)
	{ 		print SEQUENCE2 "\n$seq";			}
	else
	{  chomp($seq);
		print SEQUENCE2 "$seq";}
	}
	close(SEQUENCE1);
close(SEQUENCE2);

open (SEQUENCE1, "temp.fa") or die "error(input):$!";
open(SEQUENCE2, ">$_[0]\_length.txt") or die "error (output):$!";
# 打开sequence2.seq文件，可以向文件中写入去回车后的DNA序列。写前会删除文件以前的内容。
print SEQUENCE2 "name\tpredicted_sequence_length\n";
while(my $name = <SEQUENCE1>)
{ if($name =~ />/){
         chomp $name;
         $name =~ s/>//g;
       my $seq = <SEQUENCE1>;
             $seq =~ s/[\r\n]//g;
	my $seqlength = length($seq);
     print SEQUENCE2 "$name\t$seqlength\n";
	}
        }
close (SEQUENCE1);
close (SEQUENCE2);
unlink "temp.fa";
}




sub get_targetscan_format
{
my ($file,$ref,$taxid)=@_;
open (LIST, "$file.txt") or die "error(input1):$!";
open(TEMP, ">$file\_targetscan.txt") or die "error (output1):$!";
print "$file\_targetscan.txt\n";
open(MISS, ">$file.mis") or die "error (output2):$!";
 print TEMP "RefSeqID\tGeneID\tGeneSymbol\tSpeciesID\tUTR_sequence\n";
my %sublist = ();
my $line;
while ( $line= <LIST>)
{     chomp $line;
      my @name = split(/\t/,$line);
      my $name = $name[0];
      $sublist{$name} = 1;
}
my $count = keys(%sublist);
             # print "there are $count unique record for mirTarget!\n";

my %gotlist = ();
foreach my $source(@{$ref})
{
open (FASTA, "$source") or die "error(input2):$!";
      print "$source\n";
my $name;
my $seq;
my $fastaline;
OUTER: while ( $fastaline= <FASTA> )
{
        chomp $fastaline;
	$name="KKKKKKK";
	if($fastaline =~ />/)
	{         $fastaline =~ />([\w\:\|\_\.\-\+\d]+)/;
                  $name = $1;
                  #print "name:$name\n";
                  $seq = '';
                 if ( (exists($sublist{$name}))  and (not exists($gotlist{$name})))
                 {

	            	  while($fastaline= <FASTA>)
	           	    {
                                    if(not $fastaline =~ />/)
	               		   {
                                    $seq .= $fastaline;
                                    }
	          		   else
                                    {
                                    $seq = &fa_out($seq);
                                    print TEMP "$name\t$name\t$name\t$taxid\t$seq\n";
                                    $gotlist{$name} = 1;
                                    $seq = '';
                                    redo OUTER;
                                    }
         		    }#while
	         }#if
        }##if
}#while OUT
                          if($name ne "KKKKKKK")
                          {
                            if ( (exists($sublist{$name}))  and (not exists($gotlist{$name})))
                               {
                                  $seq = &fa_out($seq);
                                  print TEMP "$name\t$name\t$name\t$taxid\t$seq\n";
                                  $gotlist{$name} = 1;
                                 $seq = '';
                               }
                           }
close FASTA
}#foreach

$count = keys %gotlist;
    print "\n\nthere are $count circRNAs for mirTarget!\n";
foreach(keys %sublist)
{
     my $tag = $_;
       if(not exists($gotlist{$tag}))
         {
           print MISS "$tag","\n";
         }#if
}#foreach
 close (LIST);
 close (TEMP);
 close (MISS);
}


sub fasta_out
{
              my $j=0;
	      my $temp = $_[0];
              $temp =~ s/\W//g;
              my $fasta = '';
my $length = length $temp;
my $count = $length / 60 ;
my $rest = $length % 60;
for($j=0;$j<=$count;$j++)
{
my $line = substr($temp,60*$j,60);
                         if( $line )
                          {
                          $fasta .=  $line."\n";
                          }
}
              return ($fasta);
}

sub fa_out
{
              my $j=0;
	      my $temp = $_[0];
              $temp =~ s/\W//g;
              return ($temp);
}

sub anti_reverse{
              my $string=$_[0];
                 $string =~ s/\W//g;
                 $string=reverse($string);
                 $string=~ tr/ATUCGatucg/TAAGCtaagc/;
                 return($string);
}

sub circRNA_forTargetScan_out
{
              my $j=0;
	      my $temp = $_[0];
              $temp =~ s/\W//g;
              my $fasta = '';
my $length = length $temp;
    my $add='';

    if($length > 120000)
    {
       my $left = substr($temp,0,60000);
       my $right = substr($temp,$length-60000,60000);
          $temp=$left.$right;
          $add= substr($temp,0,20);
    }
    elsif($length > 20 )
    {
       $add= substr($temp,0,20);
       $length+=20;
    }
    else
    {
       $add= $temp;
       $length +=$length;
    }
     $temp=$temp.$add;
     $length = length $temp;
my $count = $length / 60 ;
my $rest = $length % 60;
for($j=0;$j<=$count;$j++)
{
my $line = substr($temp,60*$j,60);
                         if( $line )
                          {
                          $fasta .=  $line."\n";
                          }
}
              return ($fasta);
}




sub get_sub_list
{
  if($_[0] eq "--help")
  {
     print "\&get_sub_list(file1 -anh file2 -bnh -cnh -dwww -pwww)\n";

     print "    -anh\n";
     print "        n >= 0 : represents the key column of file1, 0-based.\n";
     print "        /h/    : file1 has head, otherwise file1 doesn't have head.\n\n";
     print "    -bnh is similar to a1 to config file2\n\n";
     print "    -cnh like nh\n";
     print "        n>0       : output the complement file, otherwise only output the sublist file.\n";
     print "        c1 =~ /h/ : output file should have head from file2.\n";
     print "        c1 =~ /m/ : output the losted element in file1.\n\n";
     print "    -pwww is the postfix of output file name. Default \"_contain_list\".\n\n";
     print "    -dwww is the output file name. if -d is on, -p is masked.\n";
     print "Example\n\n    \&get_sub_list(\$file1,\"-a0h\",\$file2,\"-b5h\",\"1mh\",\"-pannotation\")\n";
     print "    \&get_sub_list(\$file1,\"-a0h\",\$file2,\"-b5h\",\"-c1mh\",\"-pannotation\")\n\n";
  }
  else
  {
        my $input1="";
        my $input2="";
        my $para1="0";
        my $para2="0";
        my $para3="0";
        my $para4="_contain_list";
        my $para5="";

        my $num = @_;
        my $tag = 0;

        if($num < 2)
        {
               print "the input parameters is too little!\nPlease refer the help below:\n\n";
               &get_sub_list("--help");
               return;
        }
        else
        {
               $tag++;
        }

        my $j;
         for(my $i=0;$i < @_;$i++)
            {
               if(-e "$_[$i]")
               {
                   $input1=$_[$i];
                   $j= $i+1;
                   last;
               }
            }
         for(my $i=$j; $i < @_; $i++)
            {
            # print "$i\n";
                if(-e "$_[$i]")
               {
                   $input2=$_[$i];
                   last;
               }
            }

         unless($input1 && $input2)
           {
                print "can\'t find two exist files!\nPlease refer the help below:\n\n";
                &get_sub_list("--help");
                return;
           }

            for(my $i=0;$i<@_;$i++)
            {
               if($_[$i] =~ /\-a/)
               {
                   $para1=$_[$i];
               }
               elsif($_[$i] =~ /\-b/)
               {
                    $para2=$_[$i];
               }
               elsif($_[$i] =~ /\-c/)
               {
                    $para3=$_[$i];
               }
               elsif($_[$i] =~ /\-d/)
               {
                     $para5=$_[$i];
               }
               elsif($_[$i] =~ /\-p/)
               {
                      $para4=$_[$i];
               }
            }


        if($tag)
        {
          my $col1=0;
          my $col2=0;
          my $complement=0;

                     if( $para1 =~ /\d+/)
                     {
                       $col1 = $&;
                     }
                     if($para2 =~ /\d+/)
                     {
                       $col2 =$&;
                     }
                     if( $para3 =~ /\d+/)
                     {
                       $complement = $&;
                     }

                     my $output1;

                         if($input1 =~ /\.\w+$/)
                         {
                            $output1 =$`.$para4;
                         }

                     if($para5)
                     {
                         $output = $para5;
                     }
                     else
                     {
                         $output = $output1.".txt";
                     }
                open (TOTAL, "$input2") or die "error(input2):$!";
                open (LIST, "$input1") or die "error(input2):$!";
                open (CONTAIN, ">$output") or die "error (output1):$!";

                if($complement)
                {
                            open (UNIQ, ">$output1\_complement_list.txt") or die "error(input2):$!";
                }
                 if($para3 =~ /m/)
                {
                            open (MISS, ">$output1\_miss.txt") or die "error(input2):$!";
                }

                my %sublist = ();
                my $line;
                while ( $line= <LIST>)
                {     chomp $line;
                      my @name = split(/\t/,$line);
                      my $name = $name[$col1];
                         $sublist{$name} = 1;
                }

                my $count = keys(%sublist);
                my %gotlist = ();
                while ( $line= <TOTAL> )
                {

                  if( $line =~ /^\#/)
                   {  print "$line";   }
                   else
                   {    my $name = "KKKKKKK";
                        my @terms = split(/\t/, $line);
                        $name = $terms[$col2];
                        chomp $name;
                       	if ( exists($sublist{$name}))
                        {
                           $gotlist{$name} = 1;
                           print CONTAIN "$line";
                          }
                         else
                         {
                           if($complement)
                           {
                              print  UNIQ "$line";
                           }
                         }
                    }

                }



                my $count1 = keys %gotlist;
                print "there are $count unique records in the list file!\n";
                print "there are $count1 unique records in the gotlist file!\n";


                if($para3 =~ /m/)
                {
                    my $m=0;
                      foreach(keys %sublist)
                     {
                          my $tag = $_;
                            if(not exists($gotlist{$tag}))
                              {
                                $m++;
                                print MISS "$tag","\n";
                              }
                     }
                     print "$m records are not found in the target file!\n";
                }

                close TOTAL;
                close LIST;
                close CONTAIN;
                if($para3 =~ /m/)
                {
                 close MISS;
                }
                if($complement)
                {
                   close UNIQ;
                }
        }### if num == 2
  }##if else --help
  1;
}##sub

#########################################
sub input_table_by_poss
{
my ($database,$filename,$table)=@_;
&change_end($filename);
copy "$filename.txt.tmp","$shortdir\\$filename.txt";
open FILE1, ">input_table_by_poss.sql";
print FILE1 "
use $database
drop table $table
go
select A.name as gene_id,A.name as gene_name,A.name as transcript_id,
A.name as transcript_name,A.name as gene_biotype,A.name as source,A.name as protein_id,
A.name as uniprotID,A.name as GeneSymbol,A.name as IsoformName,A.name as proteinAcc,A.name as canonical
into $table
from UCSC_tables..poss as A
go";
close FILE1;
                   system "sqlcmd -S $server -U $user -P $pw -i input_table_by_poss.sql";
                   system("bcp $database.dbo\.$table in $shortdir/$filename.txt -c -S $server -U $user -P $pw");
unlink "input_table_by_poss.sql";
unlink "$shortdir\\$filename.txt";
unlink "$filename";
}

sub input_poss_inf_sql
{
my ($database,$filename,$table)=@_;
&change_end($filename);
copy "$filename.txt.tmp","$shortdir\\$filename.txt";
open FILE1, ">input_table_by_poss.sql";
print FILE1 "
use $_[0]
drop table $_[2]
go
select A.name as transID,A.*,A.name as GeneName,A.name as Source,A.txStart as SourceID,A.txStart as RNAlength
into $_[2]
from UCSC_tables..poss as A
go";
close FILE1;
                   system "sqlcmd -S $server -U $user -P $pw -i input_table_by_poss.sql";
                   system("bcp $database.dbo\.$table in $shortdir/$filename\.txt -c -S $server -U $user -P $pw");
unlink "input_table_by_poss.sql";
unlink "$shortdir\\$filename.txt";
unlink "$filename";
}

sub input_poss_sql
{
my ($database,$filename,$table)=@_;
&change_end($filename);
copy "$filename.txt.tmp","$shortdir\\$filename.txt";
open FILE1, ">input_table_by_poss.sql";
print FILE1 "
use $_[0]
drop table $_[2]
go
select A.*
into $_[2]
from UCSC_tables..poss as A
go";
close FILE1;
                   system "sqlcmd -S $server -U $user -P $pw -i input_table_by_poss.sql";
                   system("bcp $database.dbo\.$table in $shortdir/$filename\.txt -c -S $server -U $user -P $pw");
unlink "input_table_by_poss.sql";
unlink "$shortdir\\$filename.txt";
unlink "$filename";
}


sub get_table_from_sql
{
	my ($database,$table)=@_;
my $cmd ="bcp $database.dbo.$table out $shortdir/$table.txt -c -S $server -U $user -P $pw";
print "$cmd\n"; system $cmd;
open (INPUT, "$shortdir/$table.txt") or die "error(can't open $shortdir/$table.txt):$!";
open(OUTPUT, ">$table.txt") or die "error (can't create $table.txt):$!";
while(<INPUT>)
{
    chomp; chomp;
    print OUTPUT "$_\n";
}
close INPUT;
close OUTPUT;
unlink("$shortdir/$table.txt");
print "successful get a table!\n";
}


sub get_table_with_head
{
	my ($database,$table)=@_;
open FILE1, ">get_table_head.sql";
print FILE1 "
use $database
go
if OBJECT_ID(\'$table\_column_name\',\'U\') is not null
drop table $table\_column_name
go
select column_name
into $table\_column_name
from $database.INFORMATION_SCHEMA.COLUMNS
where TABLE_NAME=\'$table\'
go";
close FILE1;
system "sqlcmd -S $server -U $user -P $pw -i get_table_head.sql>get_table_head.log";
unlink "get_table_head.sql";
  &get_table_from_sql($database,"$table\_column_name");
  &get_table_from_sql($database,$table);
  &drop_table($database,"temp");
  rename "$table.txt","temp_$table.txt";

open (INPUT, "temp_$table.txt") or die "error(can't open temp_$table.txt):$!";
open (INPUT1, "$table\_column_name.txt") or die "error(can't open $table\_column_name.txt):$!";
open(OUTPUT, ">$table.txt") or die "error (can't create $table.txt):$!";
my @head;
while(<INPUT1>)
{
chomp;
push @head,$_;
}
my $head = join("\t",@head);
print OUTPUT "$head\n";

while(<INPUT>)
{
print OUTPUT $_;
}
close INPUT;
close INPUT1;
close OUTPUT;
unlink "temp_$table.txt";
unlink "$table\_column_name.txt";
print "successful get a table!\n";
}


sub creat_database
{
 my ($database,$force)=@_;
print "\n\n\n\t\t\tcreat_database.sql\n";
open FILE1, ">creat_database.sql";
print FILE1 "
if not EXISTS(SELECT * FROM sysdatabases WHERE name='$database')
create database $database;
";
close FILE1;
                   system "sqlcmd -S $server -U $user -P $pw -i creat_database.sql";
                   print "\n\n\n\t\t\tcreat_database.sql successful\n"        	
}

sub drop_table
{
           my ($database,$table)=@_;
           open FILE1, ">drop_table.sql";
print FILE1"
use $database
go
if OBJECT_ID(\'$table\',\'U\') is not null
	 begin
	 drop table $table
	 end
go
";
           close FILE1;
 system "sqlcmd -S $server -U $user -P $pw -i drop_table.sql";
 unlink("drop_table.sql");
}



         # &input_table_auto_head($database,$file,$table);
sub input_table_auto_head
{
my ($database,$filename,$table)=@_;
open (INPUT, "$filename.txt") or die "error(can't open $filename.txt):$!";
open(OUTPUT, ">$filename.txt.tmp") or die "error (can't create $filename.txt.tmp):$!";

print "read the head:\t";
my $line = <INPUT>;
chomp $line;
my @title = split(/\t/,$line);
my %title;
foreach (@title)
{
  s/[\+\*\#\[\]\-\(\)\{\}\?\<\>\,\;\:\"\']+/ /g;
  s/[\s\.]/_/g;
  s/_+/_/g;
  s/^_//;
  s/_$//;
  s/\//_/;
  $_='c'.$_,if(/^\d/);
  my $t = $_;
     $t =~ s/[\w\d]+//g;
     if($t)
     {
           print "title is not good $_\t$t\n";
           return;
     }
  if(lc($_) eq 'end')
  {
    $_ ="txEnd";
  }
  if(exists($title{$_}))
  {
    $_ .="_1";
    $title{$_}=1;
  }
  else
  {
    $title{$_}=1;
  }
}
print join("\t",@title),"\n";
my @type;
my @type_tag;

for(my $i=0;$i<@title;$i++)
   {
     $type_tag[$i]=1;   ###1 int 2 float 3 varchar255 4 varchar 8000
   }
print "determine the column :\n";
while($line=<INPUT>)
{
chomp $line;
my @line = split(/\t/,$line);
   for(my $i=0;$i<@title;$i++)
   {
     if($line[$i])
     {
       if(length($line[$i]) > 255)
       {
        $type_tag[$i]=4;
        }
     }
     if($line[$i] && $line[$i] =~ /\D/)
     {

       if(isnumeric($line[$i]))
       {
           $type_tag[$i]=2;
        }
       elsif($type_tag[$i] <3)
        {
          $type_tag[$i]=3;
        }
     }
   }
}
seek(INPUT,0,0);
$line=<INPUT>;
while($line=<INPUT>)
{
   chomp $line;
   my @line = split(/\t/,$line);
   for(my $i=0;$i<@title;$i++)
   {
     unless($line[$i])
     {
         if($type_tag[$i]<=2)
         {
           $line[$i]=0;
         }
         else
         {
           $line[$i]='';
         }
     }
   }
   $line = join("\t",@line[0..(@title-1)]);
   print OUTPUT "$line\n";
}


close INPUT;
close OUTPUT;

for(my $i=0;$i<@title;$i++)
   {
     if($type_tag[$i] == 1)
     {
        $type[$i] = 'int';
      }
     elsif($type_tag[$i] == 2)
     {
        $type[$i] = 'decimal(38,4)';
     }
     elsif($type_tag[$i] == 3)
     {
        $type[$i] = 'varchar(255)';
     }
     elsif($type_tag[$i] == 4)
     {
        $type[$i] = 'varchar(max)';
     }
   }
rename "$filename.txt.tmp","$shortdir\\$filename.txt";

for(my $i=0;$i<@title;$i++)
   {
     $title[$i] .= " $type[$i] null";   ###1 int 2 float 3 varchar255 4 varchar 8000
   }
my $title = join(", ", @title);

open FILE1, ">input_table_auto_head.sql";
print FILE1 "
use $database
go
if OBJECT_ID(\'$table\',\'U\') is not null
begin
drop table $table
end
go
create table $table ($title)
go
";
close FILE1;
system "sqlcmd -S $server -U $user -P $pw -i input_table_auto_head\.sql>input_table_auto_head.log";
system("bcp $database.dbo.$table in $shortdir/$filename.txt -c -S $server -U $user -P $pw");
 &check_sql_log("input_table_auto_head.log");
 #unlink "input_table_auto_head.sql";
 unlink "$shortdir/$filename.txt";
 print "successful_update a table!\n\n";
}

sub check_sql_log
{
  my $log=shift;
  open(LOG, "<$log") or die "error (can't create $log):$!";	
  while(<LOG>)
  {
	  die("there is an error in the sql log file\n"),if(/error/);
	  }
  close LOG;	
}

sub isnumeric
{
    my $val = shift;
    return length( do { no warnings "numeric"; $val & "" } ) > 0;
}

sub checkNumber{
    return shift =~ /^[+\-]?([1-9]\d*|0)(\.\d+)?([eE][+\-]?([1-9]\d*|0)(\.\d+)?)?$/;
}


             ##   &input_table_sql($database,$filename,$table,\@title,\@type);
sub input_table_sql
{
my ($database,$filename,$table,$title_ref,$type_ref)=@_;
print "\n\n\n";
print "\t\t\tInput table: $filename...............................\n";
print "\n\n\n";

&change_end($filename);
rename "$filename.txt.tmp","$shortdir/$filename.txt";
my @title;
for(my $i=0;$i<@$title_ref;$i++)
   {
     $title[$i] = "$title_ref->[$i] $type_ref->[$i] null";   ###1 int 2 float 3 varchar255 4 varchar 8000
   }
my $title = join(", ", @title);


open FILE1, ">input_table_sql.sql";
print FILE1 "
use $database
if OBJECT_ID(\'$table\',\'U\') is not null
begin
drop table $table
end
go
create table $table ($title)
go";
close FILE1;
                   system "sqlcmd -S $server -U $user -P $pw -i input_table_sql.sql>input_table_sql.log";
                   system("bcp $database.dbo\.$table in $shortdir/$filename\.txt -c -S $server -U $user -P $pw");
unlink "input_table_sql.sql";
unlink "$shortdir/$filename.txt";
print "\n\n\n";
print "\t\t\tSuccessful input table: $filename\n";
print "\n\n\n";

}

sub change_end
{

        if(-f "$_[0].txt")
        {
         open (INPUT, "$_[0].txt") or die "error(OUTPUT $_[0].txt.tmp):$!";
         open (OUTPUT, ">$_[0].txt.tmp") or die "error(OUTPUT $_[0].txt.tmp):$!";
         while(<INPUT>)
         {
            chomp;
            s/[\r\n]//g;
             print OUTPUT "$_\n";
         }
         close INPUT;
         close OUTPUT;
        }
        else
        {
               print "can't find file: $_[0].txt";
        }
}

sub check_column
{
open INPUT, "$_[0].bed";
open OUTPUT, ">$_[0]";
while(<INPUT>)
{
  chomp;
  s/[\r\n]//g;
  my @line = split(/\t/,$_);
  for(my $i=0;$i<$_[1];$i++)
  {
     if(not defined $line[$i])
          {
           $line[$i]='';
           print "$i","kkk\n";
           }
  }
  my $line = join("\t",@line);
  print OUTPUT "$line\n";
}
close INPUT;
close OUTPUT;
print "check the columns!\n";
}

sub input_bed
{
my($database,$filename,$table,$column)=@_;
my @titles=qw(chrom txStart txEnd name score strand thickStart thickEnd itemRgb blockCount blockSizes blockStarts);
my @title_use=@titles[0..$column-1];
my $title=join(",A.",@title_use);
print $title,"\n";
&check_column($filename,$column);
copy "$_[1]","$shortdir\\$_[1].txt";
open FILE1, ">input_bed_table.sql";
print FILE1 "
use $database
if OBJECT_ID('$table','U') is not null
         begin
         drop table $table
         end
go
select A.$title
into $table
from UCSC_tables..bed as A
go
";
close FILE1;
                   system "sqlcmd -S $server -U $user -P $pw -i input_bed_table.sql";
                   system("bcp $database.dbo\.$table in $shortdir/$filename\.txt -c -S $server -U $user -P $pw");
unlink "input_bed_table.sql";
unlink "$shortdir\\$filename.txt";
unlink "$filename";
}

sub psl_delete
{
open INPUT, "$_[0].psl";
open OUTPUT, ">$_[0].psll";
  my $line = <INPUT>;
  $line = <INPUT>;
  $line = <INPUT>;
  $line = <INPUT>;
  $line = <INPUT>;
print OUTPUT "$line", while($line=<INPUT>);
close INPUT;
close OUTPUT;
}

sub input_psl
{
	my($database,$filename,$table)=@_;
 &psl_delete($filename);
copy "$_[1].psll","$shortdir/$_[1].txt";
open FILE1, ">input_psl.sql";
print FILE1 "
use $database
drop table $table
select A.*
into $table
from UCSC_tables..psl as A
go";
close FILE1;
       system "sqlcmd -S $server -U $user -P $pw -i input_psl.sql";
       system("bcp $database.dbo.$table in $shortdir/$filename.txt -c -S $server -U $user -P $pw");
unlink "input_psl.sql";
unlink "$shortdir/$filename.txt";
unlink "$filename.psll";
}

sub input_blast
{
	my($database,$filename,$table)=@_;
 &change_end($filename);
 copy "$filename.txt.tmp","$shortdir\\$filename.txt";
open FILE1, ">input_blast.sql";
print FILE1 "
use $database
drop table $table
go
select A.*
into $table
from UCSC_tables..blast as A
where A.qname=\'kkk\'
go";
close FILE1;
                   system "sqlcmd -S $server -U $user -P $pw -i input_blast.sql";
                   system("bcp $database.dbo\.$table in $shortdir/$filename\.txt -c -S $server -U $user -P $pw");
 unlink "input_blast.sql";
 unlink "$shortdir\\$filename.txt";
 unlink "$filename";
}

sub input_novoalign_native
{
      my($database,$file,$table)=@_;
        &novoalign_native_treat($file);
        &input_table_auto_head($database,$file,$table);
}

sub novoalign_native_treat
{
open (INPUT, "$_[0].native") or die "error(can't open $_[1].txt):$!";
open(OUTPUT, ">$_[0].txt") or die "error (can't create temp.txt):$!";
my $line;
my @title=qw(readID readtype readseq basequalities status alignmentscore alignmentquality alignmentseq alignstart strand pairseq xxx1 xxx2);
$line =join("\t",@title);
print OUTPUT $line,"\n";
 while($line=<INPUT>)
{
   print OUTPUT $line, if($line !~ /^#/);
}
close INPUT;
close OUTPUT;
}

sub input_sam
{
    my ($database,$file,$table)=@_;
    &sam_treat($file);
    &input_table_auto_head($database,$file,$table);
}
sub sam_treat
{
open (INPUT, "$_[0].sam") or die "error(can't open $_[0].sam):$!";
open(OUTPUT, ">$_[0].txt") or die "error (can't create temp.txt):$!";
my $line;
my @title=qw(qname qname_flag tname position mappingquality CIGAR RNEXT PNEXT TLEN sequence qualities add1 add2 add3);
$line =join("\t",@title);
print OUTPUT $line,"\n";
 while($line=<INPUT>)
{
   print OUTPUT $line, if($line !~ /^\@PG/ && $line !~ /^\@SQ/ && $line !~ /^\@HD/);
}
close INPUT;
close OUTPUT;
}

sub input_gbCdnaInfo
{
	 my ($database,$filename,$table)=@_;
 &check_column($filename,22);
 copy "$filename","$shortdir\\$filename.txt";
open FILE1, ">input_gbCdnaInfo.sql";
print FILE1 "
use $database
drop table $table
go
select A.*
into $table
from UCSC_tables..gbCdnaInfo as A
where A.qname=\'kkk\'
go";
close FILE1;
                   system "sqlcmd -S $server -U $user -P $pw -i input_gbCdnaInfo.sql";
                   system("bcp $database.dbo.$table in $shortdir/$filename.txt -c -S $server -U $user -P $pw");
 unlink "input_gbCdnaInfo.sql";
 unlink "$shortdir\\$filename.txt";
 unlink "$filename";
}


sub time_interval
{
     my ($time1,$time2)=@_;
     my @time=split(/:/,$time1);
     my $s1 = timelocal($time[5],$time[4],$time[3],$time[2],$time[1],$time[0]);
        @time=split(/:/,$time2);
     my $s2 = timelocal($time[5],$time[4],$time[3],$time[2],$time[1],$time[0]);
     my $s = abs($s1 - $s2);
     if($s>=86400)
     {
         my $day = int($s/86400);
         my $hour = int(($s-$day*86400)/3600);
         my $minut = int(($s-$day*86400-$hour*3600)/60);
         return("$day:$hour:$minut");
     }
     elsif($s>=3600)
     {
         my $hour = int($s/3600);
         my $minut = int(($s-$hour*3600)/60);
         return("$hour:$minut");
     }
     else
     {
         my $minut = int($s/60);
         return("$minut");
     }
}

sub gettime {

my ($sec, $min,$hour,$mday, $mon, $year) = localtime(time);
# ($sec, $min, $hour, $mday, $mon, $year,undef,undef,undef) = localtime(time);
# $year += 1900;
# $mon  = sprintf("%02d", $mon+1);
    # printf ("%02d:%02d:%02d:%02d\n", $mday,$hour,$min,$sec);
    my $time= sprintf("%02d:%02d:%02d:%02d:%02d:%02d\n", $year,$mon,$mday,$hour,$min,$sec);
    return($time);
}
sub get_day
{
 my ($sec,$min,$hour,$mday,$mon,$year,$wday,$yday,$isdst) = localtime(time);
   $year += 1900;
   $mon++;

   if($mon<10)
   {
     $mon='0'.$mon;
    }
   if($mday<10)
   {
     $mday='0'.$mday;
   }
   my $num = $year.$mon.$mday;
   return($num);
}
sub two_excel_writer
{
	my ($sampleid,$samples,$sample_group,$fq_read_count,$circ_count,$usr,$mydir)=@_;
	my %sampleid=%{$sampleid};
	my %samples=%{$samples};
	my %fq_read_count=%{$fq_read_count};
	my %circ_count=%{$circ_count};
	my %usr=%{$usr};
	my %sample_group=%{$sample_group};
	print "$usr\t$usr->{'spn'}\n";	
	
   my $excel_new = Excel::Writer::XLSX->new('read_statistics.xlsx');
   my $worksheet = $excel_new->add_worksheet();

   $worksheet->write(0,0,"Table 1. Reads statistics");
   $worksheet->write(1,0,"Sample");
   $worksheet->write(1,1,"Raw reads");
   $worksheet->write(1,2,"Mapped reads");
   $worksheet->write(1,3,"CircRNA number");

   my $row=3;
   foreach my $sample (sort { $sampleid{$a} <=> $sampleid{$b} } keys %samples)
   {
      $worksheet->write($row-1,0,$samples{$sample});
      $worksheet->write($row-1,1,$fq_read_count{$sample}*2);
      $worksheet->write($row-1,2,$fq_read_count{$sample}*2-$unmapped{$sample});
      $worksheet->write($row-1,3,$circ_count{$sample});
      $row++;
   }
   $excel_new->close();
   ########################################################################################
    $excel_new = Excel::Writer::XLSX->new('sample & groups.xlsx');
    $worksheet = $excel_new->add_worksheet();

   $worksheet->write(0,0,"Sample Information");
   $worksheet->write(1,0,"Species");
   $worksheet->write(1,1,$usr{'spe'});
   $worksheet->write(2,0,"Sample type");
   $worksheet->write(2,1,$usr{'type'});
   $worksheet->write(3,0,"Sample number");
   $worksheet->write(3,1,$usr{'spn'});

   $worksheet->write(4,0,"Sequencing mode");
   $worksheet->write(4,1,"Paired-end");
   $worksheet->write(5,0,"Genome build");
   $worksheet->write(5,1,"HG19");
   $worksheet->write(6,0,"Sample ID");
   $worksheet->write(6,1,"Sample Name");
   $worksheet->write(6,2,"Group Name");
   $worksheet->write(6,3,"Quality Status");

   $row=1;
   foreach my $sample (sort { $sampleid{$a} <=> $sampleid{$b} } keys %samples)
   {
      $worksheet->write($row+6,0,$row);
      $worksheet->write($row+6,1,$samples{$sample});
      $worksheet->write($row+6,2,$sample_group{$samples{$sample}});
      $worksheet->write($row+6,3,"OK");
      $row++;
   }
   $excel_new->close();
}

sub two_excel_win32
{
	my ($sampleid,$samples,$sample_group,$fq_read_count,$circ_count,$usr,$mydir)=@_;
	my %sampleid=%{$sampleid};
	my %samples=%{$samples};
	my %fq_read_count=%{$fq_read_count};
	my %circ_count=%{$circ_count};
	my %usr=%{$usr};
	my %sample_group=%{$sample_group};
	print "$usr\t$usr->{'spn'}\n";
   ######################################################################   read_statistics.xlsx
my $Excel = Win32::OLE->new('Excel.Application', sub{$_[0]->Quit})
or die "Excel 初始化失败，你可能没有安装Excel！";
   $Excel->{DisplayAlerts} = 'False';    #关掉excel的提示，比如是否保存修改之类的
my $Book = $Excel->Workbooks->Add();
my $worksheet=$Book->Worksheets(1);
   $worksheet->Range("A1:D100")->Font->{Name}="Time new roman";
   $worksheet->Cells(1,1)->{Value}="Table 1. Reads statistics";
   $worksheet->Range("A1")->Font->{FontStyle}="Bold";
   $worksheet->Cells(2,1)->{Value}="Sample";
   $worksheet->Cells(2,2)->{Value}="Raw reads";
   $worksheet->Cells(2,3)->{Value}="Mapped reads";
   $worksheet->Cells(2,4)->{Value}="CircRNA number";
   $worksheet->Range("A2:D2")->Font->{FontStyle}="Bold";
   $worksheet->Range("A2:D2")->Interior->{ColorIndex}=45;
    $worksheet->Columns("A")->{ColumnWidth}=10;
    $worksheet->Columns("B:C")->{ColumnWidth}=15;
    $worksheet->Columns("D")->{ColumnWidth}=18;
   my $row=3;
   foreach my $sample (sort { $sampleid{$a} <=> $sampleid{$b} } keys %samples)
   {
      $worksheet->Cells($row,1)->{Value}=$samples{$sample};
      $worksheet->Cells($row,2)->{Value}=$fq_read_count{$sample}*2;
      $worksheet->Cells($row,3)->{Value}=$fq_read_count{$sample}*2;
      $worksheet->Cells($row,4)->{Value}=$circ_count{$sample};
      $row++;
   }

   $Book->SaveAs("$shortdir\\test.xlsx");
   $Book->close;
   
   ######################################################################   sample & groups.xlsx
my $Book1 = $Excel->Workbooks->Add();
   $worksheet=$Book1->Worksheets(1);

   $worksheet->Columns("A")->{ColumnWidth}=19;
   $worksheet->Columns("B:D")->{ColumnWidth}=18;
   $worksheet->Range("A1:D100")->Font->{Name}="Time new roman";
   $worksheet->Range("A1:D100")->Font->{Size}=11;
   $worksheet->Range("A1:A7")->Font->{FontStyle}="Bold";
   $worksheet->Range("A1")->Font->{Size}=12;
    $worksheet->Range("A2:A7")->Interior->{ColorIndex}=45;                         ### work
    $worksheet->Range("A7:D7")->Interior->{ColorIndex}=45;
    $worksheet->Range("B4")->{HorizontalAlignment}= -4131;
    $worksheet->Range("A7:D7")->Font->{FontStyle}="Bold";
   # $worksheet->Range("A1:A7")->{Interior}->{Color} = 'RGB(255,255,0)';             ### fail
   # $worksheet->Range("A1:A7")->{Interior}->{Color} = $Excel->RGB(255,255,0);     ### fail

   $worksheet->Cells(1,1)->{Value}="Sample Information";
   $worksheet->Cells(2,1)->{Value}="Species";
   $worksheet->Cells(2,2)->{Value}=$usr{'spe'};
   $worksheet->Cells(3,1)->{Value}="Sample type";
   $worksheet->Cells(3,2)->{Value}=$usr{'type'};
   $worksheet->Cells(4,1)->{Value}="Sample number";
   $worksheet->Cells(4,2)->{Value}= $usr{'spn'};
   $worksheet->Cells(5,1)->{Value}="Sequencing mode";
   $worksheet->Cells(5,2)->{Value}="Paired-end";
   $worksheet->Cells(6,1)->{Value}="Genome build";
   $worksheet->Cells(6,2)->{Value}="HG19";
   $worksheet->Cells(7,1)->{Value}="Sample ID";
   $worksheet->Cells(7,2)->{Value}="Sample Name";
   $worksheet->Cells(7,3)->{Value}="Group Name";
   $worksheet->Cells(7,4)->{Value}="Quality Status";
   $row=1;
   foreach my $sample (sort { $sampleid{$a} <=> $sampleid{$b} } keys %samples)
   {
      $worksheet->Cells($row+7,1)->{Value}=$row;
      $worksheet->Cells($row+7,2)->{Value}=$samples{$sample};
      $worksheet->Cells($row+7,3)->{Value}="$sample_group{$samples{$sample}}";
      $worksheet->Cells($row+7,4)->{Value}="OK";
      $row++;
   }
   $Book1->SaveAs("$shortdir\\test1.xlsx");
   $Book1->close;
   $Excel->{DisplayAlerts} = 'True';
   undef $Book;
   undef $Excel;
   rename "$shortdir/test.xlsx", "$mydir/read_statistics.xlsx";
   rename "$shortdir/test1.xlsx","$mydir/sample & groups.xlsx";
}



sub add_column
{
         my ($lista,$lncfather,$col1,$col2,$ref,$output,$tag)=@_;
open (LIST, "$lista.txt") or die "error($lista):$!";
open (TOTAL, "$lncfather.txt") or die "error($lncfather):$!";
open(CONTAIN, ">$output.txt") or die "error (output1):$!";
### if tag >0, the head will be added into the results.
my %sublist = ();
my $line;
if($tag)
{
         $line= <LIST>;	
         $line=~ s/[\r\n]//g;
         chomp $line;
      my @name = split(/\t/,$line);
       for(my $i=0;$i<=$#name;$i++)
         {
           unless(defined($name[$i]))
           {
                  $name[$i]='';
           }
         }
        my $name=$name[$col1];
        # print "$name\n";
      if($name)
      {
          if(defined($name))
          {
                  foreach my $key(@{$ref})
                  {
					  unless(defined($name[$key]))
                       {
                           $name[$key]='';
                       }					  
				  }
                  $sublist{$name} = join("\t",@name[@{$ref}]);
                  $sublist{'head'}=$sublist{$name};
          }
      }	
	
	}
while ( $line= <LIST>)
{     chomp $line;
	  $line=~ s/[\r\n]//g;
      my @name = split(/\t/,$line);
       for(my $i=0;$i<=$#name;$i++)
         {
           unless(defined($name[$i]))
           {
                  $name[$i]='';
           }
         }
        my $name=$name[$col1];
        # print "$name\n";
      if($name)
      {
          if(defined($name))
          {
                  foreach my $key(@{$ref})
                  {
					  unless(defined($name[$key]))
                       {
                           $name[$key]='';
                       }					  
				  }
                  $sublist{$name} = join("\t",@name[@{$ref}]);
          }
      }
}
my $count = keys(%sublist);
my %gotlist = ();
my $add_columns=@{$ref};
print "add $add_columns columns!\n";

if($tag)
{
         $line= <TOTAL>;
         chomp $line;
         $line=~ s/[\r\n]//g;
        my @terms = split(/\t/, $line);
        my $name = $terms[$col2];	
        print CONTAIN "$line\t";
        print CONTAIN "$sublist{'head'}\n";
	
}
while ( $line= <TOTAL> )
{
  if( $line =~ /^\#/)
   {  print "$line";   }
   else
   {
        chomp $line;
        $line=~ s/[\r\n]//g;
        my $name = "KKKKKKK";
        my @terms = split(/\t/, $line);
        $name = $terms[$col2];
       	if ( exists($sublist{$name}))
        {
           $gotlist{$name} = 1;
           print CONTAIN "$line\t";
           print CONTAIN "$sublist{$name}\n";
           unless($sublist{$name})
           {
             print "$name\n";
           }

        }
        else
        {
          #print CONTAIN "$line";
          print CONTAIN "$line","\t " x $add_columns,"\n";
        }
    }

}

close TOTAL;
close LIST;
close CONTAIN;
print "add columns!\n";
}

sub read_config
{
	my $config_txt=shift;
print "Reading the config file....................................................\n";
unless($config_txt)
{
my @config=<*config.txt>;
   $config_txt = $config[0];
}
unless(-f $config_txt)
{
   die("can't find config_txt file!\n");	
}
my $cfg = Config::Tiny->new;
$cfg = Config::Tiny->read($config_txt);

my($usr,$samples,$sampleid,$sample_group);

            $usr->{'name'}=$cfg->{'user_infomation'}->{'name'};
            $usr->{'prj'}=$cfg->{'user_infomation'}->{'projcet_number'};
            $usr->{'dept'}=$cfg->{'user_infomation'}->{'department'};
            $usr->{'spe'}=$cfg->{'user_infomation'}->{'species'};
            $usr->{'type'}=$cfg->{'user_infomation'}->{'sample_type'};
            $usr->{'spn'}=$cfg->{'user_infomation'}->{'sample_number'};
            # $usr{'flag'}=$cfg->{'user_infomation'}->{'XoutY_X'};
            foreach my $fq ( sort keys %{$cfg->{'sample_name'}})
            {
			    $samples->{$fq}=$cfg->{'sample_name'}->{$fq};
			}
			foreach my $id ( sort {$a <=> $b} keys %{$cfg->{'sample_group'}})
            {
			    my @terms =split("\t",$cfg->{'sample_group'}->{$id});
			    $sample_group->{$terms[0]}=$terms[1];
			    $sampleid->{$terms[0]}=int($id);
			} 
            return($usr,$samples,$sampleid,$sample_group);
}

sub txt2xlsx
{
            my $txt = shift;
            my $excel_new = Excel::Writer::XLSX->new("$txt.xlsx");
            my $worksheet = $excel_new->add_worksheet();
            open (INPUT, "$txt.txt") or die "error(input2):$!";
            my $line;
            my $columns=0;
            while ( $line= <INPUT>)
            {     chomp $line;
                  my @name = split(/\t/,$line);
                  my $column = @name;
                  $columns=$column,if($columns<$column);
            }
            seek(INPUT,0,0);
            my $row=0;
            while ( $line= <INPUT>)
            {     chomp $line;
                  my @name = split(/\t/,$line);
                  for(my $i=0;$i<=$columns;$i++)
                     {
                       unless(defined($name[$i]))
                       {
                              $name[$i]='';
                       }
                     }
                  for(my $i=0;$i<=$columns;$i++)
                  {
                     $worksheet->write($row,$i,$name[$i]);
                  }
                  $row++;
            }
            close INPUT;
            $excel_new->close();
                 #       print "The table are added into excels!\n";
}

sub txt2xlsx_win32
{
my $txt = shift;
my $app_xls = Win32::OLE->new('Excel.Application', sub{$_[0]->Quit})
or die "Excel 初始化失败，你可能没有安装Excel！";
   $app_xls->{DisplayAlerts} = 'False';
my $book = $app_xls->WorkBooks->Add;
my $booksheet= $book->Worksheets(1);
open (INPUT, "$txt.txt") or die "error(input2):$!";
my $line;
my $columns=0;
while ( $line= <INPUT>)
{     chomp $line;
      my @name = split(/\t/,$line);
      my $column = @name;
      $columns=$column,if($columns<$column);
}
seek(INPUT,0,0);
my $row=1;
while ( $line= <INPUT>)
{     chomp $line;
      my @name = split(/\t/,$line);
      for(my $i=0;$i<=$columns;$i++)
         {
           unless(defined($name[$i]))
           {
                  $name[$i]='';
           }
         }
      for(my $i=0;$i<$columns;$i++)
      {
         $booksheet->Cells($row,$i+1)->{Value}=$name[$i];
      }
      $row++;
}
close INPUT;
$book->SaveAs("$shortdir\\$txt.xlsx");
$book->close;
$app_xls->{DisplayAlerts} = 'True';
rename("$shortdir/$txt.xlsx","$mydir/$txt.xlsx");
undef $book;
undef $app_xls;  #关掉所打开的excel应用
                 #       print "The table are added into excels!\n";
}

sub get_taxid
{
  my $org = shift;
  $org = lc($org);
=pod
human          9606
mouse          10090
rat            10116
zebrafish      7955
=cut
my %taxid=();
my %build=();
$taxid{'human'}=9606;       $build{'human'}="HG19";
$taxid{'mouse'}=10090;      $build{'mouse'}="MM10";
$taxid{'rat'}=10116;        $build{'rat'}="RN5";
$taxid{'zebrafish'}=7955;   $build{'zebrafish'}="danRer10";
$taxid{'wheat'}=4565;   $build{'wheat'}="wheat_ensembl29";
          if(exists($taxid{$org}))
          {
                 print "\t\t\ttaxid:$taxid{$org}\n";
                 return($taxid{$org},$build{$org});
          }
          else
          {
                 die("unknown species, please check the org or add taxid list.");
          }

}



sub get_expression_all_txt
{
   my ($anno,$expression,$sample_count,$samples)=@_;
   open (LIST, "$anno.txt") or die "error($anno):$!";
   open (TOTAL, "$expression.txt") or die "error($expression):$!";
   open(CONTAIN, ">all.txt") or die "error (output1):$!";

my %sublist = ();
my $line;
while ( $line= <LIST>)
{     chomp $line;
      my @name = split(/\t/,$line);
      my $name = $name[3];
         $line=join("\t",@name[0..2,5..$#name]);
      $sublist{$name} = $line;
}
$sublist{'CircRNAID'} = $sublist{'name'};
my %gotlist = ();
$line= <TOTAL>;
chomp $line;
my @terms=split(/\t/, $line);
for(my $i=1;$i<=$sample_count;$i++)
{
   if(exists($samples->{$terms[$i]}))
   {
         $terms[$i]='['.$samples->{$terms[$i]}.'](raw)';
   }
   else
   {
         print "$terms[$i]\n";
   }
}
for(my $i=$sample_count+1;$i<=$sample_count*2;$i++)
{
  $terms[$i]='['.$samples->{$terms[$i]}.'](normalized)';
  # print "$terms[$i]\n";
}
$line = join("\t",@terms[0..$sample_count*2]);
print CONTAIN "$line\t$sublist{'CircRNAID'}\n";


while ( $line= <TOTAL> )
{
  if( $line =~ /^\#/)
   {  print "$line";   }
   else
   {
        my $name = "KKKKKKK";
        chomp $name;
        my @terms = split(/\t/, $line);
           $name = $terms[0];
       	if ( exists($sublist{$name}))
        {
           $gotlist{$name} = 1;
           $line=join("\t",@terms[0..$sample_count*2]);
           print CONTAIN "$line\t$sublist{$name}\n";
        }
         else
         {
              return("\n\n\nFind an un-annotated circRNA $terms[0]\n");
         }
    }

}
close TOTAL;
close LIST;
close CONTAIN;
# print "successful get the sublist!\n";
}

sub get_differentially_expressed_circRNA_list
{
&gettime;
print "get differentially expressed circRNA sequences...........................................\n";
my $file4 = shift;
print $file4,"\n";
        unless(-f "$file4")
        {
               die("Can't find file $file4!\n");
        }
my %sublist = ();
        my $parser   = Spreadsheet::ParseExcel->new();
        my $excel = $parser->parse($file4);
        for my $worksheet ( $excel->worksheets() ) {
                  my $sheetname =$worksheet->get_name();
                  print "$sheetname\n";
                  &gettime;
                  if($sheetname =~ /vs/)
                  {
                       my ( $row_min, $row_max ) = $worksheet->row_range();
                       my ( $col_min, $col_max ) = $worksheet->col_range();
                      my $rowtag;
                      my $coltag;
                      my $found=0;
                      print "excel range: $row_min\t$row_max\t$col_min\t$col_max\n";
                       for(my $i=$row_min;$i<=$row_max;$i++)
                       {
                            for(my $j=$col_min;$j<=$col_max;$j++)
                             {
                                my $cell = $worksheet->get_cell( $i, $j );
                                   next unless $cell;
                                my $value = $cell->value();
                                # print "$value\n", if($j <1);
                                if($value && $value eq 'CircRNAID')
                                {
                                   $rowtag=$i;
                                   $coltag=$j;
                                   $found++;
                                   print "CircRNAID\n";
                                   last;
                                }
                             }
                             last,if($found);
                       }
                       print "circRNAID coordinate: $rowtag\t$coltag\n";

                       for(my $i=$rowtag+1;$i<=$row_max;$i++)
                       {
                            my $cell = $worksheet->get_cell( $i, $coltag );
                                   next unless $cell;
                            my $value = $cell->value();
                               if($value)
                               {
                                       $sublist{$value}=1;
                               }
                       }
                  }## sheetname vs
              } ### sheet

my $count = keys(%sublist);
print "there are $count Differentially expressed circRNAs!\n";
my $output="Differentially expressed circRNA sequences";
open(TEMP, ">$output.txt") or die "error (output1):$!";

foreach(keys %sublist)
{
           print TEMP "$_\n";
}#foreach
close TEMP;
}


sub get_differentially_expressed_circRNA_list_win32
{
&gettime;
print "get differentially expressed circRNA sequences...........................................\n";
my $file4 = shift;
print $file4,"\n";
        unless(-f "$file4")
        {
               die("Can't find file $file4!\n");
        }
my %sublist = ();
my $app_xls = Win32::OLE->new('Excel.Application', sub{$_[0]->Quit})
or die "Excel 初始化失败，你可能没有安装Excel！";
$app_xls->{DisplayAlerts} = 'False';    #关掉excel的提示，比如是否保存修改之类的


           copy "$file4","$shortdir\\$file4";
        my $dst4_name = "$shortdir\\$file4";
        my $book4 = $app_xls->WorkBooks->Open($dst4_name);
        my $sheetcnt = $book4->Worksheets->Count();
        foreach (1..$sheetcnt)
              {
                  my $book4sheet=$book4->Worksheets($_);
                  my $sheetname = $book4sheet->{Name};
                  print "$sheetname\n";
                  &gettime;
                  if($sheetname =~ /vs/)
                  {
                      # my $LastRow = $book4sheet->UsedRange->Find({What=>"*", SearchDirection=>xlPrevious,SearchOrder=>xlByRows})->{Row};
                      my $LastRow = $book4sheet->{UsedRange}->{Rows}->{Count}; 
                      # my $LastCol = $book4sheet->UsedRange->Find({What=>"*", SearchDirection=>xlPrevious,SearchOrder=>xlByColumns})->{Column};
                      my $LastCol       = $book4sheet->{UsedRange}->{Columns}->{Count}; 
                      my $rowtag;
                      my $coltag;

                       for(my $i=1;$i<=$LastRow;$i++)
                       {
                            my $value = $book4sheet->Cells($i,1)->{Value};
                               if($value && $value eq 'CircRNAID')
                               {
                                   $rowtag=$i;
                                   last;
                               }
                       }
                             for(my $j=1;$j<=$LastCol;$j++)
                             {
                               my $value = $book4sheet->Cells($rowtag,$j)->{Value};
                               if($value && $value eq 'CircRNAID')
                               {
                                   $coltag=$j;
                                   last;
                               }
                             }
                       print "circRNA coordinate: $rowtag\t$coltag\t$LastCol\n";

                       for(my $i=$rowtag+1;$i<=$LastRow;$i++)
                       {
                            my $value = $book4sheet->Cells($i,$coltag)->{Value};
                           # print "$value\n";
                               if($value)
                               {
                                       $sublist{$value}=1;
                               }
                       }
                  }## sheetname vs
              } ### sheet
              $book4->close;
        undef $book4;
        $app_xls->{DisplayAlerts} = 'True';
        undef $app_xls;  #关掉所打开的excel应用
my $count = keys(%sublist);
print "there are $count Differentially expressed circRNAs!\n";
my $output="Differentially expressed circRNA sequences";
open(TEMP, ">$output.txt") or die "error (output1):$!";

foreach(keys %sublist)
{
           print TEMP "$_\n";
}#foreach
close TEMP;
}
sub sqlserver_table_exists
{
         my ($database,$table)=@_;
         open FILE1, ">sql_inf_test.sql";
         print FILE1 "
         use $database
         if OBJECT_ID('$table','U') is not null
begin
select distinct count(*)
from $table as A
end
         ";
         close FILE1;
         my $cmd="sqlcmd -S $server -U $user -P $pw -i sql_inf_test.sql";
         print "$cmd\n";
         my $inf=`$cmd`;
         my @line=split(/[\r\n]/,$inf);
         my $line_count=@line;
         
         if($line_count>1)
         {
             my $line = $line[5];
                $line=~ /\d+/;
             my $count =$&;
             return(1); 
	     }
	     else
	     {
			return(0); 
			 }
}

sub sqlserver_table_line_count
{
         my ($database,$table)=@_;
         open FILE1, ">sql_inf_test.sql";
         print FILE1 "
         use $database
         if OBJECT_ID('$table','U') is not null
begin
select distinct count(*)
from $table as A
end
         ";
         close FILE1;
         my $cmd="sqlcmd -S $server -U $user -P $pw -i sql_inf_test.sql";
         print "$cmd\n";
         my $inf=`$cmd`;
         my @line=split(/[\r\n]/,$inf);
         my $line_count=@line;
         
         if($line_count>1)
         {
             my $line = $line[3];
                $line=~ /\d+/;
             my $count =$&;
             return($count); 
	     }
	     else
	     {
			return(0); 
			 }
}


sub get_md5_files
{
       my ($projectid,$samples) =@_;
     open (LIST, "md5.fastq.txt") or die "error(md5.fastq):$!";
     my %md5_hash=();
     my $line;
     while($line =<LIST>)
     {
              chomp $line;
           my ( $file,$hash)=split(/\t/,$line);
                $md5_hash{$file}=$hash;
     }
     close LIST;
             open (OUTPUT, ">$projectid.md5") or die "error($projectid.md5):$!";
     foreach my $sample(sort keys %{$samples})
     {
           if(exists($md5_hash{$sample."_R1.fastq.gz"}))
           {

                print OUTPUT $md5_hash{"$sample\_R1.fastq.gz"},"\t$samples->{$sample}\_R1.fastq.gz\n";

           }
           else
           {
                          die("no md5 value for $sample\_R1.fastq.gz\n");
            }
           if(exists($md5_hash{$sample."_R2.fastq.gz"}))
           {
                         print OUTPUT $md5_hash{"$sample\_R2.fastq.gz"},"\t$samples->{$sample}\_R2.fastq.gz\n";
           }
           else
           {
                          die("no md5 value for $sample\_R2.fastq.gz\n");
            }

     }
     close OUTPUT;
}

1;
