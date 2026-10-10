#!/usr/bin/env perl
use warnings;
use strict;
use Getopt::Long;
use File::Basename;

########################################################################################
  # Descript: A simply submission task script based on batch submission system.
  # Author: Xin Du; Email: lzudux2014@gmail.com; Department: Institute of Zoology, CAS
  # Date: 2021/9/14
  # Version: v1.0
  # Modification: 2022/09/14
  # Version: v2.0
########################################################################################

my ($threads,$node,$partition,$node_num,$queue,$jobprefix,$memory,$lines,$submit,$ve,$exit,$system,$log_path,$help);

GetOptions(
    'threads|t=n'  => \$threads,
	'node|n=s' => \$node,
	'partiotion|p=s' => \$partition,
    'queue|q=s'  => \$queue,
	'node_num|nn=s' => \$node_num,
    'jobprefix|jp=s' => \$jobprefix,
    'memory|m=s' => \$memory,
    'lines|l=n'  => \$lines,
	'log_path|lp=s' => \$log_path,
    'submit|s=s'  => \$submit,
    'virtua_environment|v', => \$ve,
	'system|ss=s' => \$system,
    'exit|e=s'   => \$exit,
    'help|h' => \$help
    );
my $command = shift;
chomp(my $work_path = `pwd`);
my $help_out = "
perl $0 [options] run.sh
###       Necessary      ###
The command list file

###        Options       ###
--system|-ss    string    # SEG|PBS|SLURM|LSF|auto, default=auto.
--threads|-t   number    # Number of threads, default is 1.
--nodes|-n     string    # node name, default: system allocation.
--node_num|-nn int       # node number, default: 1 node for per task.
--queue|-q     string    # All can be used queue in your server, default: not specified.
--partition|-p string    # Specifying a computation partition, default: not specified.
--lines|-l     number    # Number of line of per splited file, default is off.
--memory|m     string    # Set memery size (option unit: M or G), default is 2G.
--jobprefix|jp string    # The prefix of job, default is work. No more than 5 characters are recommended for mutil lines task. Starting with a number is not allowed.
--submit|-s    NA|sh_dir # Directly submit all jobs(yes or no) or given ready dir name, default is no.
--virtua_environment|-v  # To submit the virtual environment(such as conda or python virtua environment and etc).
--exit|-e      string    # Whether exit the script (yes:default or no).
--log_path|-lp path      # log path, default: work path.
--help|-h                # Print help information.

###    example usage    ###
easy_submit.pl -t 30 -m 40G -s yes run.sh

###   when you have prepared submit_jobs_dir  ###
easy_submit.pl -s submit_jobs_dir -e no|yes
";

die $help_out if($help);
die $help_out if(!$command and !$submit);
die $help_out if(!$command  and !-d $submit);
my $submission_type = &auto_detect_submission_system_type;

# prepare parameters
$node_num = 1 if(!$node_num);
$threads = 1 if(!$threads);
$submit = "no" if(!$submit);
#print STDOUT "$submit\n";
$memory = "2G" if(!$memory);
$memory = uc($memory);
$jobprefix = "work" if(!$jobprefix);
if(length $jobprefix > 10){
	my @tmp = split //,$jobprefix;
	if($lines){
		$jobprefix = join "",@tmp[0..4];
	}else{
		$jobprefix = join "",@tmp[0..7];
	}
}

$exit = "yes" if(!$exit);
$log_path = $work_path if(!$log_path);
my $task_dir;
$task_dir = &Run_Core($command,$submission_type) if(!-d $submit);

if(-d $submit){
    my $all_jobs = &submit_jobs($submit, $submission_type);
	die "Task submission failed!\n" if(!$all_jobs);
    &check_state($all_jobs, $submission_type) if($exit eq "no");
    exit 0;
}elsif($task_dir and $submit eq "yes"){
	#print STDOUT "test done!\n";
    my $all_jobs = &submit_jobs($task_dir, $submission_type);
	die "Task submission failed!\n" if(!$all_jobs);
    &check_state($all_jobs, $submission_type) if($exit eq "no");
}

sub Run_Core{
    my $file = shift;
	my $sys = shift;
    my $sh = basename($file);
    $sh =~ s/\.sh$//;
    my ($sec,$min,$hour,$mday,$mon,$year,$wday,$yday,$isdst) = localtime(time);
    $year += 1900;
    $mon += 1;
    my $datetime = sprintf ("%d%02d%02d%02d%02d%02d", $year,$mon,$mday,$hour,$min,$sec);
	$sh = "STDIN" if($sh eq "-");
    my $out_dir = "$log_path/$sh\_submit_jobs.$datetime";
    my @exists = <$log_path/$sh\_submit_jobs*>;
    if(@exists >= 1){
        `rm -rf $_` for(@exists);
    }
    mkdir $out_dir;
    &build_task_scripts($sys,$file,$out_dir);
    return $out_dir;
}

sub build_task_scripts{
	my ($sys,$shell_file,$out_dir) = @_;
	if($lines){
		open(F,$shell_file) || die "Can't open the $shell_file.\n";
		my $count =0;
		my $lable = 1;
		my %save;
		while(<F>){
			chomp;
			next if(/^#|^\s*$/);
			$count++;
			if($count <= $lines){
				push(@{$save{$lable}},$_);
			}else{
				$lable++;
				push(@{$save{$lable}},$_);
				$count=1;
			}
		}
		close F;
		for my $n(sort {$a<=>$b} keys %save){
			my $ord = sprintf("%03d",$n);
			my $name = $jobprefix.$ord;
			open(O,">$out_dir/$name.$sys");
			print O &scripts_info($sys,$out_dir,$ord);
			if(@{$save{$n}} == 1){
				print O "\n$save{$n}[0] && echo \"This work is completed!\" && exit 1\necho \"This work is not completed!\"\n";
			}else{
				my $n = @{$save{$n}} - 2;
				print O join "\n",@{$save{$n}}[0..$n];
				print O join "\n${$save{$n}}[-1] && echo \"This work is completed!\" && exit 1\necho \"This work is not completed!\"\n";
			}
			close O;
		}
	}else{
		my $name = $jobprefix;
		open(O,">$out_dir/$name.$sys");
		print O &scripts_info($sys,$out_dir);
		my @line;
		for(`cat $shell_file`){
			chomp;
			next if /^#|^\s*$/;
			push(@line,$_);
		}
		if(@line >= 2){
			print O join "\n",@line[0..$#line-1];
			print O "\n$line[-1] && echo \"This work is completed!\" && exit 1\necho \"This work is not completed!\"\n";
		}else{
			print O "$line[0] && echo \"This work is completed!\" && exit 1\necho \"This work is not completed!\"\n";
		}
		close O;
	}
}

sub scripts_info{
	my ($sys,$out_dir,$num) = @_;
	my $return = "";
	my $name;
	$name = $jobprefix;
	$name = $jobprefix.$num if($num);
	if($sys eq "PBS"){
		$return = "#!/bin/bash
#PBS -N $name
#PBS -l nodes=$node_num:ppn=$threads
#PBS -l mem=$memory
#PBS -S /bin/bash
#PBS -o $out_dir/$name.$sys.out
#PBS -e $out_dir/$name.$sys.err
";
		$return .= "#PBS -V\n" if($ve);
		$return .= "#PBS -q $queue\n" if($queue);
		$return .= "#PBS -l nodes=$node\n" if($node);
		$return .= "cd $work_path\n";
	}elsif($sys eq "SGE"){
		$return = "#!/bin/bash
#\$ -cwd
#\$ -N $name
#\$ -pe smp $threads
#\$ -l vf=$memory
#\$ -S /bin/bash
#\$ -o $out_dir/$name.$sys.out
#\$ -e $out_dir/$name.$sys.err
";
		$return .= "#\$ -V\n" if($ve);
		$return .= "#\$ -q $queue\n" if($queue);
		$return .= "#\$ -l hostname=$node\n" if($node);
	}elsif($sys eq "slurm"){
		$return = "#!/bin/bash
#SBATCH --job-name $name
#SBATCH -c $threads --mem $memory
#SBATCH --get-user-env
#SBATCH -o $out_dir/$name.$sys.out
#SBATCH -e $out_dir/$name.$sys.err
#SBATCH -D $work_path
#SBATCH -N $node_num
";
		$return .= "#SBATCH --partition=$partition\n" if($partition);
		#$return .= "\n" if($ve);
		$return .= "#SBATCH -w $node\n" if($node);
	}elsif($sys eq "LSF"){
		$return = "
#!/bin/bash
#BSUB -J $name
#BSUB –n $threads
#BSUB -M $memory
#BSUB -R span[hosts=$node_num]
#BSUB -o $out_dir/$name.$sys.out
#BSUB -e $out_dir/$name.$sys.err
";
		$return .= "#BSUB -q $queue\n"if($queue);
		$return .= "cd $work_path\n";
	}else{
		die "Unknown submission system type!\n";
	}
	return $return;
}

sub submit_jobs{
    my $command_dir = shift;
	my $sys = shift;
    my %jobs_id;
    my @command_file = <$command_dir/*$sys>;
    for(@command_file){
		if($sys eq "SGE"){
	        my $job_return = `qsub $_ 2>/dev/null`;
    	    if($job_return =~ /Your\s+job\s+(\d+)/){
        	    $jobs_id{$1}++;
            	print STDOUT "The job number for $_: $1\n";
	        }else{
    	        die "Error: Check $_\n";
        	}
		}elsif($sys eq "PBS"){
			my $job_return = `qsub $_ 2>/dev/null`;
			if($job_return =~ /(\S+)/){
				$jobs_id{$1}++;
				print STDOUT "The job number for $_: $1\n";
			}else{
				die "Error: Check $_\n";
			}
		}elsif($sys eq "slurm"){
			my $job_return = `sbatch $_ 2>/dev/null`;
			if($job_return =~ /Submitted\s+batch\s+job\s+(\d+)/){
				die "Error: Check the $_\n" if(!$1);
				$jobs_id{$1}++;
				print STDOUT "The job number for $_: $1\n";
			}else{
				die "Error: Check $_\n";
			}
		}else{
			my $job_return = `bsub < $_ 2>/dev/null`;
			if($job_return =~ /(\d+)/){
				$jobs_id{$1}++;
				print STDOUT "The job number for $_: $1\n";
			}else{
				die "Error: Check $_\n";
			}
		}
    }
    return \%jobs_id;
}

sub check_state{
    my $jobs = shift;
	my $sys = shift;
    my %jobs = %$jobs;
    #print join "\n",(keys %jobs),"\n";
    my %tmp;
    my $j = 0;
    DO:while(1){
        for my $job(sort keys %jobs){
            #print "The job number: ".$job."\n";
            #`qstat -j $job`;
			if($sys eq "SGE"){
            	my $state_return = `qstat -j $job 2>/dev/null`;
	            $tmp{$job}++ if(!$state_return and !exists $tmp{$job});
			}elsif($sys eq "PBS"){
				my $state_return = `qstat -f $job 2>/dev/null`;
				$tmp{$job}++ if(!$state_return and !exists $tmp{$job});
			}elsif($sys eq "slurm"){
				my $state_return = `squeue -j $job 2>/dev/null`;
				my $line_num = $state_return =~ s/\n//g;
				$tmp{$job}++ if((!$state_return or $line_num == 1) and !exists $tmp{$job});
			}else{
				my $state_return = `bjobs -o stat -noheader \"$job\" 2>/dev/null`;
			}
        }
        last DO if(scalar(keys %tmp) == scalar(keys %jobs));
        sleep 2;
    }
}

sub auto_detect_submission_system_type{
	# in order judge
	# Slurm: if sbatch is present, assume slurm works.
	# SGE:   if qstat -a can't be performed, assume SGE works.
	# LSF:   if LSF_ENVDIR exists in the environment, assume LSF works.
	# PBS:   if qstat -a can be performed, assume PBS works.
	my $type;
	chomp(my $slurm=`which sbatch 2>/dev/null`);
	chomp(my $qsub = `which qsub 2>/dev/null`);
	chomp(my $lsf = `echo \$LSF_ENVDIR`);
	if($qsub){
		my $pbs = `qstat -a 2>/dev/null`;
		if($pbs){
			print STDOUT "Your submission system is PBS\n";
			$type = "PBS";
		}else{
			print STDOUT "Your submission system is SGE\n";
			$type = "SGE";
		}
	}elsif($lsf){
		$type = "LSF";
		print STDOUT "Your submission system is LSF\n";
	}elsif($slurm){
		$type = "slurm";
		print STDOUT "Your submission system is slurm\n";
	}else{
		print STDERR "Unknown submission system!\n";
		exit 0;
	}
	return $type;
}
