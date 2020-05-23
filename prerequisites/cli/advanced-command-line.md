Advanced Command-Line
=======================

The sed command
----------------

Let's take a look at the 'sed' command. sed (short for stream editor) is a command that allows you to manipulate character data in various ways. One useful thing it can do is substitution. First, make a directory called "advanced" to work in, for this document. If you have access to /share/workshop, then make the "advanced" directory under your username there. If you don't, just put the "advanced" directory in your home. The "$USER" variable contains your username.

    cd /share/workshop/$USER/  # or, if this fails, just 'cd ~' or 'cd'
    mkdir advanced
    cd advanced/

Let's copy over a simple file to work on:

    cp /usr/share/common-licenses/BSD .

Take a look at the file:

    cat BSD

Now, let's change every occurence of the word "Redistribution" into "Mangling":

    cat BSD | sed 's/Redistribution/Mangling/gi'

Let's break down the argument to sed (within the single quotes)... The "s" means "substitute", the word between the 1st and 2nd forward slashes (i.e. /) is the word the substitute for, the word between the 2nd and 3rd slashes is the word to substitute with, and finally the "gi" at the end are flags for global substitution (i.e. substituting along an entire line instead of just the first occurence on a line), and for case insenstivity (i.e. it will ignore the case of the letters when doing the substitution).

Note that this **doesn't** change the file itself, it is simply piping the output of the cat command to sed and outputting to the screen. If you wanted to change the file itself, you could use the "-i" option to sed:

    cat BSD
    sed -i 's/Redistribution/Mangling/gi' BSD

Now if you look at the file, the lines have changed.

    cat BSD

Another useful use of sed is for capturing certain lines from a file. You can select certain lines from a file:

    sed '4q;d' BSD

This will just select the 4th line from the file.

**CHALLENGE:**
See if you can find a way to use sed to remove all the spaces from the BSD file.

More pipes
-----------

Now, let's delve into pipes a little more. Pipes are a very powerful way to look at and manipulate complex data using a series of simple programs. First take a look at the contents of the "/home" directory:

    ls /home

These are all the home directories on the system. Now let's say we wanted to find out how many directory names begin with each letter. First we cut out the first letter of the directories:

    ls /home | cut -c1

In order to do the counting, we first need to sort the data and then send it to the "uniq" command to keep only the unique occurences of a letter. The "-c" option counts up the number of occurences:

    ls /home | cut -c1 | sort | uniq -c


Now let's look at some fastq files. Link a few files into the advanced directory:

    ln -s /share/biocore-archive/Leveau_J_UCD/RNASeq_Arabidopsis_2016/00-RawData/*/*.fastq.gz .
    
Since the files are gzipped files we need to use "zcat" to look at them. zcat is just like cat except for gzipped files:

    zcat C61_S67_L006_R1_001.fastq.gz | head

Notice that each header line has the barcode for that read at the end of the line. Let's count the number of each barcode. In order to do that we need to just capture the header lines from this file. We can use "sed" to do that:

    zcat C61_S67_L006_R1_001.fastq.gz | sed -n '1~4p' | head

By default sed prints every line. In this case we are giving the "-n" option to sed which will **not** print every line. Instead, we are giving it the argument "1\~4p", which means to print the first line, then skip 4 lines and print again, and then continue to do that.

Now that we have a way to get just the headers, we need to isolate the part of the header that is the barcode. There are multiple ways to do this... we will use the cut command:

    zcat C61_S67_L006_R1_001.fastq.gz | sed -n '1~4p' | cut -d: -f10 | head

So we are using the "-d" option to cut with ":" as the argument to that option, meaning that we will be using the delimiter ":" to split the input. Then we use the "-f" option with argument "10", meaning that we want the 10th field after the split. In this case, that is the barcode.

Finally, as before, we need to sort the data and then use "uniq -c" to count. Then put it all together and run it on the entire dataset (This will take about a minute to run):

    zcat C61_S67_L006_R1_001.fastq.gz | sed -n '1~4p' | cut -d: -f10 | sort | uniq -c

Now you have a list of how many reads were categorized into each barcode. Here is a [sed tutorial](https://www.digitalocean.com/community/tutorials/the-basics-of-using-the-sed-stream-editor-to-manipulate-text-in-linux) for more exercises.

One final thing to know is that if a program does not take input from STDIN (which is needed to use it in a pipe), but instead wants a filename, you can use a single dash by itself in place of the filename and the shell will interpret that to be input from STDIN. So it would look something like this:

<div class="output">cat FILENAME | COMMAND -f - -otheroptions | ....
</div>

**CHALLENGE:**
Find the distribution of the first 5 bases of all the reads in C61_S67_L006_R1_001.fastq.gz. I.e., count the number of times the first 5 bases of every read occurs across all reads.

Process substitution
---------------------

Next, we will cover process substitution. Process substitution is a way of using the output of some software as the input file to another software without having to create intermediate files. We will use a quality-based trimmer called "sickle". We want to do quality-based read trimming on one of our fastq.gz files, but we need to give sickle an uncompressed file as input. In order to do that, we use the "gunzip" command with the "-c" option. This unzips the file and sends the output to STDOUT, instead of unzipping the file in place which is the default (This will take a few minutes to run):

    module load sickle
    sickle se -f <(gunzip -c C61_S67_L006_R1_001.fastq.gz) -t sanger -o trimmed.fa

So we are putting the gunzip command inside parentheses with a less-than symbol like so: <(COMMAND). When we do this, the output of the COMMAND gets manipulated by the shell so that sickle thinks it is a file. Sickle then uses this "file" as the input file. Take a look at the output file:

    less trimmed.fa

Loops
------

Loops are useful for quickly telling the shell to perform one operation after another, in series. For example:

    for i in {1..21}; do echo $i >> a; done  # put multiple lines of code on one line, each line terminated by ';'
    cat a
    # <1 through 21 on separate lines>

The general form is:

<div class="output">for name in {list}; do
    commands
done
</div>

The list can be a sequence of numbers or letters, or a group of files specified with wildcard characters:

    for i in {3,2,1,liftoff}; do echo $i; done  # needs more excitement!
    for i in {3,2,1,"liftoff!"}; do echo $i; done  # exclamation point will confuse the shell unless quoted

 A "while" loop is more convenient than a "for" loop ... if you don't readily know how many iterations of the loop you want:

<div class="output">while {condition}; do
    commands
done
</div>

Now, let's do some bioinformatics-y things with loops and pipes. First, let's write a command to get the nucleotide count of the first 10,000 reads in a file. Use zcat and sed to get only the read lines of a file, and then only take the first 10,000:

    zcat C61_S67_L006_R1_001.fastq.gz | sed -n '2~4p' | head -10000 | less

Use grep's "-o" option to get each nucleotide on a separate line (take a look at the man page for grep to understand how this works):

    zcat C61_S67_L006_R1_001.fastq.gz | sed -n '2~4p' | head -10000 | grep -o . | less

Finally, use sort and uniq to get the counts:

    zcat C61_S67_L006_R1_001.fastq.gz | sed -n '2~4p' | head -10000 | grep -o . | sort | uniq -c

<div class="output"> 264012 A
 243434 C
 215045 G
    278 N
 277231 T
</div>

And, voila, we have the per nucleotide count for these reads!

We just did this for one file, but what if we wanted to do it for all of our files? We certainly don't want to type the command by hand dozens of times. So we'll use a while loop. You can pipe a command into a while loop and it will iterate through each line of the input. First, get a listing of all your files:

    ls -1 *.fastq.gz

Pipe that into a while loop and read in the lines into a variable called "x". We use "$x" to get the value of the variable in that iteration of the loop:

    ls -1 *.fastq.gz | while read x; do echo $x is being processed...; done

Add the command we created above into the loop, placing $x where the filename would be and semi-colons inbetween commands:

    ls -1 *.fastq.gz | while read x; do echo $x is being processed...; zcat $x | sed -n '2~4p' | head -10000 | grep -o . | sort | uniq -c; done

When this runs it will print the name of every single file being processed and the nucleotide count for the reads from those files.

Now, let's say you wanted to write the output of each command to a separate file. We would redirect the output to a filename, but we need to create a different file name for each command and we want the file name to reflect its contents, i.e. the output file name should be based on the input file name. So we use "parameter expansion", which is fancy way of saying substitution:

    ls -1 *.fastq.gz | while read x; do echo $x is being processed...; zcat $x | sed -n '2~4p' | head -10000 | grep -o . | sort | uniq -c > ${x%.fastq.gz}.nucl_count.txt; done

This will put the output of the counting command into a file whose name is the prefix of the input file plus ".nucl_count.txt". It will do this for every input file.


Bash Scripts
---------------

A script is a set of commands that are written into a file. This file can then be "run" as a program and it will simply execute all of the commands in it. Open a new text file using the text editor "nano":

    nano get_nucl_counts.sh

Copy and Paste the following into the file:

<div class="script">#!/bin/bash

zcat $1 | sed -n '2~4p' | head -$2 | grep -o . | sort | uniq -c
</div>

Save the file and exit. Change the permissions on the file to make it executable:

    chmod a+x get_nucl_counts.sh

Now, we can run this script giving it different arguments every time. The first argument (i.e. the first text after the script name when it is run) will get put into the variable "$1". The second argument (delimited by spaces) will get put into "$2". In this case, "$1" is the file name, and "$2" is the number of reads we want to count. So, then we can run the script over and over again using different values and the command will run based on those values:

    ./get_nucl_counts.sh I593_S85_L006_R1_001.fastq.gz 1000
    ./get_nucl_counts.sh I593_S85_L006_R1_001.fastq.gz 10000
    ./get_nucl_counts.sh C64_S70_L006_R2_001.fastq.gz 555
    ./get_nucl_counts.sh <(gzip -c BSD) 10

We can also put loops into a script. We'll take the loop we created earlier and put it into a file, breaking it up for readability and using backslashes for line continuation:

    nano get_nucl_counts_loop.sh

Put this in the file and save it:

<div class="script">#!/bin/bash

ls -1 *.fastq.gz | \
while read x; do \
    echo $x is being processed...; \
    zcat $x | sed -n '2~4p' | head -$1 | \
        grep -o . | sort | uniq -c > ${x%.fastq.gz}.nucl_count.txt; \
done
</div>


Make it executable:

    chmod a+x get_nucl_counts_loop.sh

And now we can execute the entire loop using the script. Note that there is only one argument now, the number of reads to use:

    ./get_nucl_counts_loop.sh 100

Find
-----

Find is a very powerful command that is used to recursively find files/directories in a file system. First take a look at the find man page:

    man find

Notice there are LOTS of options. The simplest version of the find command will simply list every file and directory within a path, as far down as it can go:

    find /share/biocore/joshi/projects/genomes

If we want to refine the command to only show you files that end in ".fa" (i.e. fasta files), we use the "-name" option:

    find /share/biocore/joshi/projects/genomes -name "*.fa"

One of the most powerful uses of find is to execute commands on every file it finds. To do this, you use the "-exec" option. When you use that option, everything after the "-exec" is assumed to be a command, and you use the "{}" characters to substitute for the file names that it finds. So in the command below, "wc -l" will get executed sequentially for every file it finds. Finally, the exec option needs to end with a semi-colon, however, since the semi-colon is a special character that the shell will try to interpret, you need to "escape" the semi-colon with a backslash, to indicate to the shell that the semi-colon is NOT to be interpreted and just sent as is to the find command:

    find /share/biocore/joshi/projects/genomes -name "*.fa" -exec wc -l {} \;

You will probably want to Ctrl-C out of this because it will take a long time to go through them all.

**HARD CHALLENGE:**
Use the "find" command to find all the files ending in ".pm" in the /software/perl-libs directory. Then, pipe that to a while loop to grep for all occurences (case insensitive) of the word "blast" in those files. The grep should also output the file name. You will probably have to look at the man page for grep.

Xargs
------

xargs is another command that can be very useful for running a program on a long list of files. For example, the "find" commands we ran above could be run using xargs like this:

    find /share/biocore/joshi/projects/genomes -name "*.fa" | xargs wc -l

This is taking the output of the find command and then creating a list of all the filenames which it adds to the command given to xargs. So, in this case, after "xargs" comes "wc -l"... so "wc -l" will get run on the entire list of filenames from find.


.bashrc/.bash_profile, aliases & the PATH variable
-----------------------------------------------------

On a Linux system, there is usually a user-modifiable file of commands that gets run every time you log in. This is used to set up your environment the way that you want it. On our systems, the file is ".bash_profile" and it resides in your home directory. Sometimes the file is called ".bashrc" as well. Take a look at a .bash_profile:

    cat /share/workshop/.bash_profile

This one has a lot of stuff in it to set up the environment. Now take a look at your own .bash_profile. One of the things you set up was adding to the PATH variable. One thing that is very useful to add to the PATH variable is the "." directory. This allows you to execute things that are in your current directory, specified by ".". So, let's use nano to edit our .bash_profile and add "." to PATH:

    nano ~/.bash_profile

Add ":." to the PATH variable by adding this line:

**export PATH=$PATH:.**

 Now, next time you log in, "." will be in your PATH.

Another thing that is very useful are aliases. An alias is a user-defined command that is a shortcut for another command. For example, let's say you typed the command "ls -ltrh" a lot and it would be easier to have it be a simpler command. Use an alias:

    alias lt='ls -ltrh'

Now, you've created an alias that lists the contents of a directory with extra information (-l), in reverse time order of last modified (-r and -t), and with human readable file sizes (-h). Try it out:

    lt

Typing alias by itself give you a list of all the aliases:

    alias

You can now put the alias command for lt in your .bash_profile and you will have it automatically when you log in.

More grep
----------

Grep is a very powerful tool that has many applications. Grep can be used to find lines in a file that match a pattern, but also you can get lines above and below the matching line as well. The "-A" option to grep is used to specify number of lines after a match and the "-B" option is for number of lines before a match. So, for example, if you wanted to find a particular sequence in a fastq file and also the 3 other lines that form that entire fastq record, you would do this:

    zcat C61_S67_L006_R1_001.fastq.gz | grep -B1 -A2 CACAATGTTTCTGCTGCCTGAACC

This looks for the sequence "CACAATGTTTCTGCTGCCTGAACC" in the fastq file and then also prints the line before and two lines after each match. 

Another thing grep can do is regular expressions. Regular expressions are a way of specifying a search pattern. Two very useful characters in regular expressions are "^" and "$". The "^" symbol specifies the beginning of a line and the "$" specifies the end of a line. So, for example, if you wanted to find just the lines that began with "TTCCAACACA" you would do this:

    zcat C61_S67_L006_R1_001.fastq.gz | grep ^TTCCAACACA

Without the "^", grep will find any line that has "TTCCAACACA" *anywhere* in the line, not just the beginning. Conversely, if you wanted to find the lines that ended in "TAAACTTA":

    zcat C61_S67_L006_R1_001.fastq.gz | grep TAAACTTA$

There are also extended regular expression that grep can use to do more complex matches, using the "-E" option:

    zcat C61_S67_L006_R1_001.fastq.gz | grep -E '^TTCCAACACA|TAAACTTA$'

This command will find any line that begins with "TTCCAACACA" **OR** ends with "TAAACTTA". The "\|" character means OR.

**CHALLENGE:**
Find a way to use grep to match any line that has between 7 and 16 'A's in a row at the end of the line. You will probably need to look at the man page for grep.

The Prompt
------------

The Prompt is the part of the command line that is the information on the screen before where you type commands. Typically, it has your username, the name of the machine you are on, and your current directory. This, like everything else in Linux, is highly customizable. Your current prompt probably has your username, your hostname, and your current directory. However, there are many other things you could add to it if you so desired. The environment variable used for your prompt is "PS1". So to change your prompt you just need to reassign PS1. There are some special escape characters that you use to specify your username, hostname, etc... these can be found in the "PROMPTING" section of the bash man page:

    man bash

Then type "/" and "PROMPTING" to find the section. You'll see that for username it is "\u", for hostname it is "\h", and for the full current working directory it is "\w". So you would change your prompt by doing this:

    export PS1="\u@\h:\w\\$ "

The convention is to put the "@" symbol and the ":" symbol to delineate the different parts, but you can use anything you want. Also, it is convention to put a "$" at the end to indicate the end of the prompt.

You can also do all kinds of fancy things in your prompt, like color and highlighting. Here is an example of such a prompt:

    export PS1='\[\033[7;29m\]\u@\h\[\033[0m\]:\[\e[1m\]\w\[\e[m\]$ '

When you have a prompt you like, you can put it in your .bash_profile/.bashrc so that it is automatically set when you log in.

Nohup
------

The nohup (short for "no hangup") command is useful for running a job from a terminal and then wanting to exit the terminal. When you run a job, even if you put it in the background (i.e. by using "&"), the job is tied to the terminal you are running on. When you log out of that terminal, any job tied to that terminal will be killed. This is not desirable, so you can use the nohup command which will disconnect the proccess from the terminal. Simply put "nohup" in front of the command, and you will probably want to add the "&" at the end to put it in the background so you can get your prompt back. It would look something like this:

<div class="output">nohup YOUR COMMAND &
</div>

Awk
----

Awk is a simple programming language that can be used to do more complex filtering of data. Awk has many capabilities, and we are only going to touch on one of them here. One really useful thing is to filter lines of a file based on the value in a column. Let's get a file with some data:

    wget https://raw.githubusercontent.com/ucdavis-bioinformatics-training/2019_August_UCD_mRNAseq_Workshop/master/cli/DMR.GBM2.vs.NB1.bed

Take a look at the beginning of the file:

    head DMR.GBM2.vs.NB1.bed

Let's say we wanted to get only the lines where the pvalue (column 10) was below a certain value. In awk, the default delimiter is the tab character, which most bioinformatics files use as their delimiter. Using the tab as a delimiter, it assigns each value in a column to the variables $1, $2, $3, etc... for as many columns as there are. So if we wanted to get lines where the pvalue column was under 0.00000005 pvalue:

    cat DMR.GBM2.vs.NB1.bed | awk '$10 < 0.00000005'

And lines where the pvalue >= 0.00000005:

    cat DMR.GBM2.vs.NB1.bed | awk '$10 >= 0.00000005'

You can also use it to extract lines with a particular word in a column:

    cat DMR.GBM2.vs.NB1.bed | awk '$1 == "chr3"'

A double equals (==) is used for equality comparisons. This will pull out lines where the chromosome column is "chr3".

Take a look at the [awk manual](https://www.gnu.org/software/gawk/manual/gawk.html) to learn more about the capabilities of awk.

**HARD CHALLENGE**:
Go through the list of genomes (as in the Find section) and this time only search down a maximum of 6 directories and also follow symbolic links in the search. Then extract only those files that are part of either the zebrafish or C. elegans genomes. For each of those files, get the number of characters in the file and then only print files whose character count is less than 10000. You will have to probably use find, grep, xargs, wc, and awk. You will need to look at the manual pages for each of those commands. You should be able to do this just using pipes and the commands (i.e. no intermediate files).

