#!/share/apps/bin/perl

### Import packages used by the programm.
use strict;
use warnings;
use Getopt::Long;
use DBD::mysql;
use Data::Dumper;
use File::Basename;
use vars qw($help $output $out_log);
use Bio::SeqIO;
use Bio::Tools::GFF;
use Bio::SeqFeature::Generic;
use Bio::Location::Split;
use Bio::SearchIO;

### Define program's options.
my $VERSION = '2.1';
my $lastmodif = '2014/06/28';
my $time = localtime;
my $userName = 'root';
my $password = '********';
my $hostName = 'localhost';
my $port = 3306;
my $dbName = '3b';
my $output;
my $outLog;

&GetOptions ( "h|help"           => \&help,
              "u|username=s"     => \$userName,
              "p|password=s"     => \$password,
              "host=s"           => \$hostName,
              "d|dbname=s"       => \$dbName,
              "o"                => \$output,
              "l"                => \$outLog );

@ARGV or &help;
&main(\@ARGV);

### This program work by passing through different subroutines, this is the main sub.
sub main {
	my $self = {};
	bless $self;
	$self->{inputFiles} = shift;
	$self->setOptions();
	$self->connectToDB();
	$self->startLog() if ($self->{outLog});
	foreach (@{$self->{inputFiles}}){
		$self->{currentFile} = $_;
		unless (-e $self->{currentFile}) {
			print STDERR "*** Could not find " . $self->{currentFile} . " ***\n";
			next;
		}
		$self->readFile();
	}
	$self->{dbi}->disconnect();
	exit(0);
}

### This sub allow the programm to assign variables as options.
sub setOptions {
	my $self = shift;
	$self->{userName} = $userName;
	$self->{hostName} = $hostName;
	$self->{password} = $password;
	$self->{port} = $port;
	$self->{dbName} = $dbName;
	$self->{output} = $output;
	$self->{outLog} = $outLog;
	$self->{log}->{time} = $time;
	return 1;
}

### Sub used to connect mySQL database (containing all TE infos).
sub connectToDB {
	my $self = shift;
	$self->{dbi} = DBI->connect("DBI:mysql:database=" . $self->{dbName} . ";host=" . $self->{hostName}, $self->{userName}, $self->{password})
	                            or die("Could not connect to database" . "\n" . DBI->errstr);
	return 1;
}

### Put log's file values at zero.
sub startLog {
	my $self = shift;
	$self->{log}->{noTE} = 0 ;
	$self->{log}->{"oneTE"} = 0 ;
	$self->{log}->{"multiTE"} = 0 ;
	$self->{log}->{"TEtot"} = 0 ;
	$self->{log}->{"smallTE"} = 0 ;
	$self->{log}->{"junctionTOT"} = 0 ;
	$self->{log}->{"junctionOK"} = 0 ;
	$self->{log}->{"junctionN"} = 0 ;
	$self->{log}->{"junctionOV"} = 0 ;
	$self->{log}->{"st_spNeg"} = 0 ;
	return 1;
}

### Open file containing reference FASTA sequences.
sub readFile {
	my $self = shift;
	$self->{seqIoObj} = Bio::SeqIO->new(-format => "fasta", 
	                                    -file => $self->{currentFile});
	if($self->{output}){
		$self->{outSeqIoObj} = Bio::SeqIO->new(-format => "fasta", 
		                                       -file => ">" . $self->{currentFile} . ".junctions.fna");
	}
	else{
		$self->{outSeqIoObj} = Bio::SeqIO->new(-format => "fasta", 
		                                       -fh => \*STDOUT);
	}
	$self->{log}->{countScaff}=0;
	while ( $self->{seqObj} = $self->{seqIoObj}->next_seq() ) {
		$self->{log}->{countScaff}++;
		$self->{scaffId} = $self->{seqObj}->display_id();
	}
		print STDERR "# " . $self->{scaffId} . "\n";
		
		### mySQL request.
		$self->{sql} = 	"SELECT id," . 
		                       "TEid," .
		                       "start," .
		                       "stop," . 
		                       "start-75 AS st5," . 
		                       "start+75 AS sp5," .
		                       "stop-75 AS st3," .
		                       "stop+75 AS sp3," .
		                       "scaff" .
		                 "FROM annotTE3b" .   
		                 "WHERE scaff='" . $self->{scaffId} . "'" .
		                 "ORDER BY start" ;

		### Execute mySQL request.
		$self->{sth} = $self->{dbi}->prepare($self->{sql});
		$self->{sth}->execute();
		undef $self->{junction};
		while(my $row = $self->{sth}->fetchrow_hashref){
			if ($row->{stop}-$row->{start}+1 < 75){
				$self->{log}->{"smallTE"}++;
				next;
			}
			$self->{junction}->{$self->{scaffId}}->{$row->{id}.'_5'}->{idTe}  = $row->{id}; 
			$self->{junction}->{$self->{scaffId}}->{$row->{id}.'_5'}->{scaff} = $row->{scaff}; 
			$self->{junction}->{$self->{scaffId}}->{$row->{id}.'_5'}->{start} = $row->{st5}; 
			$self->{junction}->{$self->{scaffId}}->{$row->{id}.'_5'}->{stop}  = $row->{sp5}; 

			$self->{junction}->{$self->{scaffId}}->{$row->{id}.'_3'}->{idTe}  = $row->{id}; 
			$self->{junction}->{$self->{scaffId}}->{$row->{id}.'_3'}->{scaff} = $row->{scaff}; 
			$self->{junction}->{$self->{scaffId}}->{$row->{id}.'_3'}->{start} = $row->{st3}; 
			$self->{junction}->{$self->{scaffId}}->{$row->{id}.'_3'}->{stop}  = $row->{sp3}; 
		}

		### Total junctions number in scaffold analysed.
		$self->{log}->{"TEtot"} = scalar(keys(%{$self->{junction}->{$self->{scaffId}}}));
		
		### Count TE number per scaffold (0, 1, several)
		if    ($self->{log}->{TEtot} == 0) { $self->{log}->{"noTE"}++; }
		elsif ($self->{log}->{TEtot} == 1) { $self->{log}->{"oneTE"}++; }
		else                               { $self->{log}->{"multiTE"}++; }
		
		### Define $junction as scaff's key
		print STDERR "->filterOverlappingJunctions\n";
		$self->filterOverlappingJunctions();
		print STDERR "->getSeq\n";
		$self->getSeqOfJunctions();
		print STDERR "->printLog\n";
		$self->printLog() if($self->{outLog});
	}
} 

sub filterOverlappingJunctions {
	my $self = shift;
	my $kRef = 0; 
	REF:foreach my $idRef (sort {$self->{junction}->{$self->{scaffId}}->{$a}->{start} <=> 
	                             $self->{junction}->{$self->{scaffId}}->{$b}->{start}} 
						   keys(%{$self->{junction}->{$self->{scaffId}}})){
		### In case this junction has been discarded based on overlap.
		if($self->{result}->{$self->{scaffId}}->{discard}->{$idRef}) {
			next;
		}
		if ($self->{junction}->{$self->{scaffId}}->{$idRef}->{start} <= 0){
			$self->{result}->{$self->{scaffId}}->{discard}->{$idRef}='startOfScaff';
			next;
		}
		elsif ($self->{junction}->{$self->{scaffId}}->{$idRef}->{stop} > $self->{seqObj}->length){
			$self->{result}->{$self->{scaffId}}->{discard}->{$idRef}='endOfScaff';
			next;
		}

		### Declare refID as "kept" junction before checking for overlapping junctions.
		### Use {start} as value to be able to sort on start in &getSeqOfJunction.
		$self->{result}->{$self->{scaffId}}->{keep}->{$idRef} = $self->{junction}->{$self->{scaffId}}->{$idRef}->{start};
		$kRef++;
		my $startRef = $self->{junction}->{$self->{scaffId}}->{$idRef}->{start};
		my $stopRef = $self->{junction}->{$self->{scaffId}}->{$idRef}->{stop};
		my $kComp = 0;
		COMP:foreach my $idComp (sort {$self->{junction}->{$self->{scaffId}}->{$a}->{start} <=> 
		                               $self->{junction}->{$self->{scaffId}}->{$b}->{start}} 
								keys(%{$self->{junction}->{$self->{scaffId}}})){ 
			my $startComp = $self->{junction}->{$self->{scaffId}}->{$idComp}->{start};
			my $stopComp = $self->{junction}->{$self->{scaffId}}->{$idComp}->{stop};
			
			### Next if same junction ID.
			if ($idRef eq $idComp){next;} 
			
			### Last the foreach loop if start_comp > stop_ref (i.e juction are overlapping).
			if ($startComp > $stopRef){last;}
			
			### Test the overlap.
			if ($startRef<=$stopComp and $startRef>=$startComp){
				$self->{result}->{$self->{scaffId}}->{overlap}->{$idRef}->{++$kComp} = $idComp;
				$self->{result}->{$self->{scaffId}}->{discard}->{$idComp}='overlap';
			}
			elsif ($stopRef<=$stopComp and $stopRef>=$startComp){
				$self->{result}->{$self->{scaffId}}->{overlap}->{$idRef}->{++$kComp}=$idComp;
				$self->{result}->{$self->{scaffId}}->{discard}->{$idComp}='overlap';
			}
		}
	}

	if(scalar(keys%{$self->{result}->{$self->{scaffId}}->{keep}}) + 
       scalar(keys%{$self->{result}->{$self->{scaffId}}->{discard}}) != 
       scalar(keys%{$self->{junction}->{$self->{scaffId}}}) ) {
		die "Error: junctions keep+discard differ from totalNumberOfJunction\n";
	}
	return 1;
}

### Get sequences of junction from fasta files and coordonates from database.
sub getSeqOfJunctions {
	my $self = shift;
	### In $rslt hash, non-overlapping sequences are put inside $subseq.
	foreach my $id ( sort {$self->{result}->{$self->{scaffId}}->{keep}->{$a} <=> 
	                       $self->{result}->{$self->{scaffId}}->{keep}->{$b} }
	                 keys%{$self->{result}->{$self->{scaffId}}->{keep}}){
		my $subseq = $self->{seqObj}->subseq($self->{junction}->{$self->{scaffId}}->{$id}->{start},
		                                     $self->{junction}->{$self->{scaffId}}->{$id}->{stop});
		### No N content test.
		if ($subseq =~/N/i){ 
			$self->{log}->{junctionN}++;
			$self->{result}->{$self->{scaffId}}->{withNs}->{$id}=1;
			next;
		}
		else {
			my $subseqObj = Bio::Seq ->new('-seq' => $subseq, 
			                               '-id' => $id,
			                               '-description' => $self->{scaffId} . "_" . 
			                                                 $self->{junction}->{$self->{scaffId}}->{$id}->{start} . "," .
			                                                 $self->{junction}->{$self->{scaffId}}->{$id}->{stop}); 
			                                
			$self->{log}->{junctionOK}++;
			$self->{result}->{$self->{scaffId}}->{withoutNs}->{$id}=1;
			$self->{outSeqIoObj}->write_seq($subseqObj); 
		}
	}
	return 1;
}

### Print final log to resume the analysis.
sub printLog {
	my $self = shift;

	print join("\t", $self->{scaffId}, "totalNbTE", $self->{log}->{TEtot}), "\n";
	print join("\t", $self->{scaffId}, "totalNbSmallTE", $self->{log}->{smallTE}), "\n";
	print join("\t", $self->{scaffId}, "nbDiscard", scalar keys%{$self->{result}->{$self->{scaffId}}->{discard}}),"\n";
	print join("\t", $self->{scaffId}, "nbKeep", scalar keys%{$self->{result}->{$self->{scaffId}}->{keep}}), "\n";
	print join("\t", $self->{scaffId},  "nbWithNs", scalar keys%{$self->{result}->{$self->{scaffId}}->{withNs}}), "\n";
	print join("\t", $self->{scaffId}, "nbWithoutNs", scalar keys%{$self->{result}->{$self->{scaffId}}->{withoutNs}});

	foreach my $id (sort keys%{$self->{result}->{$self->{scaffId}}->{discard}}){
		print join("\t", $self->{scaffId}, "discardId", $id, $self->{result}->{$self->{scaffId}}->{discard}->{$id}), "\n";
	}
	return 1;
}

### Program's help.
sub help {
my $prog = basename($0);
print STDERR <<EOF ;

# # # $prog # # #
#
# CREATED:    2014/04/22
# LAST MODIF: $lastmodif
# AUTHOR:     Emeric ******** (INRA Clermont-Ferrand)
# VERSION:    $VERSION
# PURPOSE:    This script is used to design non overlapping junctions 
#             between TE/TE or TE/LowCopyDNA
#             with 75 nt around junctions from annotated mySQL database.
# USAGE:
#             $prog  [OPTIONS]  fastaFile1  fastaFile2  ...
#
# OPTIONS:
#             -h       print this help
#             -l       print details in log file (XLS format)
#             -o       redirect output into a file [default: STDOUT]
#             -u       <string>  sql username [default: root]
#             -p       <string>  sql password [default: ********]
#             -host    <string>  sql hostname [default: localhost]
#             -d       <string>  sql DB name  [default: 3b]

EOF
	exit(1) ;
}
