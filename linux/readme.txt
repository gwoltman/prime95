
Welcome to the Great Internet Mersenne Prime Search!

To use this program you must agree to the prize rules at
http://www.mersenne.org/prize.htm

In case you ever forget, the URL is http://www.mersenne.org/prime.htm.
My email address is woltman@alum.mit.edu.  For networking questions,
contact Scott Kurowski's technical support team at primenet@mersenne.org.


FILE LIST
---------

readme.txt	This file.
mprime		The program to factor and run Lucas-Lehmer tests on
		Mersenne numbers.
whatsnew.txt	A list of new features in mprime.
stress.txt	A discussion of issues relating to stress testing a computer.
undoc.txt	A list of formerly undocumented and unsupported features.
prime.ini	A file containing your preferences.  The menu choices
		and dialog boxes are used to change your preferences.
local.ini	Like prime.ini, this file contains more preferences.
		The reason there are two files is discussed later.
worktodo.ini	A list of exponents the program will be factoring
		and/or Lucas-Lehmer testing.
results.txt	Mprime writes its results to this file.
prime.log	A text file listing all messages that have been sent
		to the PrimeNet server.
prime.spl	A binary file of messages that have not yet been sent to
		the PrimeNet server.
pNNNNNNN &	Intermediate files produced by mprime to resume
qNNNNNNN	computation where it left off.
mNNNNNNN	Intermediate files produced during P-1 factoring.


WHAT IS THIS PROGRAM?
---------------------

This program is used to find Mersenne Prime numbers.  See
http://www.utm.edu/research/primes/mersenne.shtml for a good
description of Mersenne primes.  Mersenne numbers can be proved
composite (not prime) by either finding a factor or by running
a Lucas-Lehmer primality test.


COMMAND LINE ARGUMENTS
----------------------

These command line arguments can be used to schedule mprime to
run only at certain times of the day or at a different priority.  Note
that raising the program's priority will not make it run any faster
on an idle machine.

-An		This is used to run two or more copies of mprime
		from the same directory.  Using this command line argument
		causes mprime to use a different set of filenames for the
		INI files, the results file, the log file, and the spool file.
		Just use a different value of n for each extra copy of
		mprime you start.
-Bn		Run in background.  This is the same as adding an "&" at
		the end of the mprime command line.  The optional value n
		is the number of CPUs on your machine.  One mprime will
		be run in background for each CPU.  "mprime -B2" is equivalent
		to "mprime -B" followed by "mprime -B -A1".
-c		Contact the PrimeNet server then exit.  Useful, for
		scheduling server communication as a cron job or
		as part of a script that dials an ISP.
-d		Prints more detailed information to stdout.  Normally
		mprime does not send any output to stdout.
-m		Bring up the menus to set mprime's preferences.
-t		Run the torture test.  Same as Options/Torture Test.
-v		Print the version number of mprime.
-Wdirectory	This tells mprime to find all its files in a different
		directory than the executable.


INSTRUCTIONS
------------

There are two ways to use this program.  The automatic way uses
a central server, which we call the PrimeNet server, to get work to do
and report your results.  Anyone with Internet access, even dial-up users,
should use this method.  You do not need a permanent connection to the
Internet.

The second method is the manual method.  It requires a little more work
and monitoring.  I recommend this for computers with no Internet access
or with some kind of firewall problem that prevents the automatic method
from working.

If you are running this program at your place of employment, you should
first GET PERMISSION from your network administrator or boss.  This is
especially true if you are installing the software on several machines.
Some companies are reluctant to run any software they are not familiar with.


INSTRUCTIONS FOR THE AUTOMATIC METHOD
-------------------------------------

1)  Download and decompress mprime.tgz.  You've probably done this already
    since you are reading this file.
2)  Connect to the Internet.
3)  Run mprime.  You will be see 3 sets of questions:
3a) The first set gathers information about you.  Enter your name and email
    address.  Optionally enter a user ID, password, and computer ID.  If you
    are using several computers, use the same user ID and password but a
    unique computer ID on each machine.  An easy-to-remember user ID will be
    helpful if you plan to visit the PrimeNet server's web page to view
    reports on your progress.  If you do not enter a user ID or if you pick a
    user ID that is already in use, then the server will assign a user ID.
3b) The second set of questions gathers information about your machine.
    Fill in roughly how many hours a day you leave your computer running.
3c) The third set of questions gathers information about using the
    Primenet server.  Note that mprime will never dial-up to connect to
    the Internet, rather it waits for a time when you are already connected
    to contact the server.  Mprime will now contact the PrimeNet server to
    get some work for your computer to do.
4)  If you are using a proxy server you may need to do some further
    configuration.  See the later secion on "SETTING UP A PROXY SERVER".
5)  Configure your startup scripts so that mprime runs every time you
    boot your computer.


MANUAL METHOD INSTRUCTIONS
--------------------------

1)  Use the Web (http://mersenne.org/ips/manualtests.html) to create
    a userid for yourself and to get a set of exponents to work on.
    Alternatively, you can get some work from
    http://www.mersenne.org/range2.htm.  Copy these exponents to a
    file called worktodo.ini
2)  Run mprime.  You will be see 3 sets of questions:
2a) The first set gathers information about you.  Enter your name and email
    address.
2b) The second set of questions gathers information about your machine.
    Fill in roughly how many hours a day you leave your computer running.
2c) The third set of questions gathers information about using the
    Primenet server.  Answer "No" to the first question.
3)  Configure your startup scripts so that mprime runs every time you
    boot your computer.
4)  Once a month or when done with your exponents, use the web pages
    again to send the file "results.txt" to the PrimeNet server.
    It is important to do this so the exponents you are testing are
    not reassigned to someone else.


NOTES
-----

Let mprime run at all times.  It runs at low priority, making use of all
your idle CPU cycles.  It should not interfere with your normal work.
Let the program run overnight and on weekends.  Never turn
your computer off.  Turn off your monitor to conserve energy.  NOTE:
Running your computer non-stop could increase your electric bill by $30
per year or more.

It can take many CPU days to test a large Mersenne number.  This program
can be safely interrupted by using ^C or kill or shutdown to write
intermediate results to disk.  This program also saves intermediate
results to disk every 30 minutes in case there is a power failure.

To fully utilize a dual-CPU machine, you must run two copies of
mprime.  Run one copy of mprime as described above.  Run the second
copy of mprime with the -A1 switch.  Make sure your startup scripts
start both executables.

You can compare your computer's speed with other users by checking the
site http://www2.tripnet.se/~nlg/mersenne/benchmk.htm.

You can get a report of your PrimeNet server account status on the web
(http://mersenne.org/primenet/status.shtml).  Type your UserID and
password into the web form, and click "Get Account Report".  You can
get your UserID and password from the Test/User Information
dialog box.

Before testing an exponent, the program may perform an hour long
self-test to make sure the Lucas-Lehmer code is running properly
on your system.  If this hour long test is interrupted it will restart
from the beginning when the program resumes testing.

If you have overclocked your machine, I recommend running the torture
test for a couple of days.  The longer you run the torture test
the greater the chance that you will uncover an error caused by
overheating or overstressed memory.

Depending on the exponent being tested, the program may decide that it
would be wise to invest some time checking for small factors before
running a Lucas-Lehmer test.  Furthermore, the program may start factoring
exponents before a previous Lucas-Lehmer test completes.  This is normal!
The program will resume the Lucas-Lehmer test when the factoring
completes.

You can configure this program to have different properties at different
times of the day and/or to not run during certain times of the day.
Unfortunately, you must manual edit the prime.ini file.  Let's say you
want to install the program on a friend's machine and he runs a screen
saver at night.  He also runs a disk defragmenter at midnight on weekdays.
This prime.ini file will run the program at a higher priority than his
screen saver at night and on weekends.  It also sleeps for an hour when
his defragmenter starts running.  Finally, at night fewer save files are
generated to allow his disk to stay powered down longer.
	UserID=foo
	Time=1-5/8:30-17:30
	Priority=1
	DiskWriteTime=30
	Time=1-5/1:00-8:30,1-5/17:30-24:00,6-7/0:00-24:00
	Priority=5
	DiskWriteTime=240
The 1-5 refers to days of the week, with Monday=1 and Sunday=7.  The 
time portion refers to the hours of the day based on a 24-hour clock.
You do not need to specify days of the week (e.g. Time=1-7/0:00-8:00
is the same as Time=0:00-8:00).  Unpredictable results will occur if
there are overlapping time intervals.  Also note that any options that
appear in the Time= sections should not appear earlier in the prime.ini
file and you can no longer edit these options from theuser interface.


SETTING UP A PROXY SERVER
-------------------------

For the latest information on this and other networking issues visit
the FAQ at http://mersenne.org/ips/faq.html

Create a file called "primenet.ini" in the same directory as mprime.
Add the line "ProxyHost=foo.com:1234" where foo.com is the proxy and
1234 is the port (which defaults to 8080 if not specified). This line must
be in the .ini file ONLY if connection via a proxy is actually required.


SETTING AVAILABLE MEMORY
------------------------

The P-1 factoring step prior to running a Lucas-Lehmer test is more
effective if it is given more memory to work with.  However, if you let
the program use too much memory then the performance of ALL programs will
suffer.  The good news is that 98% of the time the program uses less
than 8MB.  In fact, the program will work just fine if you instruct the
program to use only 8MB or less.

So how do you intelligently choose the available memory settings?  Below
are some steps you might take to figure this out:

1)  Be conservative.  It is better to set the available memory too low
than too high.  Setting the value too high can cause thrashing which
slows down all programs.  Remember, the program will only use the
extra memory in stage 2 of P-1 factoring (about 12 hours a month).

2)  Start with how much memory is installed in your machine.  Allow a
reasonable amount of memory for the OS and whatever background tasks
you run (say 16 or 24MB).  This represents the maximum value you should use.
The program won't let you enter more than 90% of installed memory.

3)  Assuming you run your machine 24 hours a day, what hours of the
day do you not use your computer?  Make these your nighttime hours and
let the program use a lot of memory during these hours.  But reduce this
value if you also run batch jobs at night.

4)  Factor in the information below about minimum, reasonable, and
desirable memory amounts for some sample exponents.

	Exponent	Minimum		Reasonable	Desirable
	--------	-------		----------	---------
	 6000000	 12MB		   23MB		  33MB
	10000000	 19MB		   36MB		  53MB
	33000000	 65MB		  125MB		 185MB

For example, my machine is a dual-processor with 128MB of memory.
Each CPU runs an LL test on an exponent near 10 million.  On the off
chance that both CPUs wind up doing P-1 factoring at the same time,
I don't want to set the available memory to more than half of the
128MB.  I guess Windows NT can survive on 24MB of memory.  Thus, I set
the available memory to (128 - 24) / 2 or 52MB.  This is my nighttime setting.
During the day, I set the available memory to 24MB.  I can always stop
prime95 if it is doing P-1 factoring and I detect memory thrashing.  More
casual users will probably want to set the daytime memory to 8MB so they
don't have to worry about it.

If at all in doubt, leave the settings at 8MB.  The worst that will
happen is you end up running a Lucas-Lehmer primality test when stage 2
of P-1 factoring would have found a factor.


PROGRAM OUTPUT
--------------

On screen you will see:

Factoring M400037 to 2^54 is 3.02% complete. Clocks: 24235224=0.121 sec.
	This means mprime is trying to find a small factor of 2^400037-1.
	It is 3.02% of the way though looking at factors below 2^54.  When
	this completes it may start looking for factors less than 2^55.
Iteration: 941400 / 1667747 [56.45%].  Per iteration time: 0.109 sec. (21889762 clocks)
	This means mprime just finished the 941400th iteration of a
	Lucas-Lehmer primality test.  The program must execute 1667747
	iterations to complete the primality test.  The average iteration
	took 21889762 "clock cycles" or 0.109 seconds.

The results file and screen will include lines that look like:

M2645701 has a factor: 13412891051374103
	This means to 2^2645701-1 is not prime.  It is divisible
	by 13412891051374103.
M2123027 no factor to 2^57, WV2: 14780E25
	This means 2^2123027-1 has no factors less than 2^57.  The Mersenne
	number may or may not be prime.  A Lucas-Lehmer test is needed
	to determine the primality of the Mersenne number.  WV2 is
	the program version number.  14780E25 is a checksum to guard
	against email transmission errors.
M1992031 is not prime. Res64: 6549369F4962ADE0. WV2: B253EF24,1414032,00000000
	This means 2^1992031-1 is not prime - a Lucas-Lehmer test says so.
	The last 64 bits of the last number in the Lucas-Lehmer sequence
	is 6549369F4962ADE0.  At some future date, another person will verify
	this 64-bit result by rerunning the Lucas-Lehmer test.  WV2 is the
	program	version number.  B253EF24 is a checksum to guard against email
	transmission errors.  1414032 can be ignored it is used as part
	of the double-checking process.  The final 00000000 value is a set
	of 4 counters.  These count the number of errors that occurred during
	the Lucas-Lehmer test.
M11213 is prime! WV2: 579A579A
	This means 2^11213-1 is a Mersenne prime!  WV2 is the program
	version number.  579A579A is a checksum to guard against email
	transmission errors.


RUNNING MPRIME ON SEVERAL COMPUTERS
-----------------------------------

The easiest way to do this is to first set up mprime on one computer.
Next copy all the files to the second computer.  Delete the local.ini
file and worktodo.ini files.  These files contain information that
is specific to the first computer.  Start mprime on the second
computer and optionally use Test/User Information to
give the second computer a unique computer ID.  Repeat this process
for all the computers you wish to run mprime on.

If you do not follow the instruction above, be sure you use 
Test/User Information to give each computer the same
userid and password.  Failure to do this will result in all 
your work being "credited" to different user IDs.


TEST MENU
---------

The PrimeNet menu choice lets you configure the type of work you get
from the server.  The "Use PrimeNet..." option can be turned on to
switch from the manual method to the automatic method.
The "Request whatever..." box should be left checked.  However, if you
are running a slow computer and don't mind waiting several months for
a single Lucas-Lehmer test to complete OR you are running a faster
computer and would rather do factoring, then uncheck this box and
choose a different type of work to do.  The "Always have at least
this many days of work queued up" value should be changed based on
how often you connect to the Internet.  As long as you connect at least
once in the given time period, mprime will have an uninterrupted stream
of work.  However, the program will not checkout more than 20 exponents
no matter what this value is set at.

The User Information menu choice lets you change your name and email
address.  Your name will be used in credits and "top producers"
web pages.  Your email address may be used to send you a
newsletter once every few months.  Select the "Receive occasional
newsletters..." option to enable these emails.  The userid you
use will appear in the status pages on the PrimeNet server.  You can
also set the ComputerID field if you have several computers and want
to keep track of which computers produced which results.  You can create
a team by choosing the "Create a team..." option.  The user name you
entered becomes the team name and a message is sent to the server which
prevents the team name from changing.  Thus, you can handout the userid
and password for others to join your team and these team members are not
able to accidentally or maliciously change the team name that appears
on the statistics page.

The Vacation/Holiday menu choice lets mprime update the expected completion
dates on the PrimeNet server.  This will prevent one of your exponents
from being reassigned while you're gone.  Also, if you are taking a
long vacation, mprime will get extra exponents to test.  For example,
if you are leaving immediately for a 6 week vacation, connect to the
Internet and choose 42 days.

The Status menu choice will tell you what exponents you are working on.
It will also estimate how long that will take and your chances of finding
a new Mersenne prime.

The Continue menu choice lets you resume mprime after you have stopped it.

The Exit menu choice lets you exit the program.


ADVANCED MENU
-------------

You should not need to use the Advanced menu.  This menu choice is
provided only for those who are curious to play with.

The Test choice can be used to run a Lucas-Lehmer test on one Mersenne
number.  Enter the Mersenne number's exponent - this must be a prime
number between 5 and 79300000.

The Time choice can be used to see how long each iteration of a Lucas-Lehmer
test will take on your computer and how long it will take to test a
given exponent.  For example, if you want to know how long a Lucas-Lehmer
test will take to test the exponent 876543, choose Advanced/Time and
enter 876543 for 100 iterations.  A Pentium will automatically time
and compute the results.  486 owners will need to time how long the operation
takes by hand, divide that time by 100 and multiply by 876543 to compute how
long a full Lucas-Lehmer test will take.

The ECM choice lets you factor small Mersenne exponents using the
Elliptic Curve Method of factoring.  Select a few exponents and bounds
to factor from the http://www.mersenne.org/ecm.htm web page.
Note:  You do not reserve exponents to work on, several people can
do ECM factoring on the same exponent.  The program uses a random
number generator to select elliptic curves to test.  You must email
results to me at woltman@alum.mit.edu - primenet does not support ECM
factoring.  Thus, ECM factoring of F20 which uses a 64K FFT will use
a minimum of 192 * 64K or 12MB of memory.  You can also edit the
worktodo.ini file directly.  For example:
	ECM=751,3000000,0,100,0,0,0,0
The first value is the exponent.  The second value is bound #1.  The
third value is bound #2 - leave it as zero.  The fourth value is the
number of curves to test.  The fifth value is no longer used. 
The sixth value is the specific curve to test - it is only used in
debugging.  The seventh value is 0 for 2^N-1 factoring, 1 for 2^N+1
factoring.  The eighth value is no longer used.

The P-1 choice lets you factor Mersenne numbers using the P-1
method of factoring.  The web pages at http://www.mersenne.org/ecm.htm
will have more information about coordinating the P-1 factoring effort.
When P-1 completes it produces a save file that others can use to do further
P-1 factoring.  You must email results to me at woltman@alum.mit.edu -
primenet does not support P-1 factoring.  You can also edit the worktodo.ini
file directly.  For example:
	Pminus1=751,3000000,0,0,0
The first value is the exponent.  The second value is bound #1.  The
third value is bound #2 or zero for default.  The fourth value is 0 for
2^N-1 factoring, 1 for 2^N+1 factoring.  The fifth value is no longer used.

The Priority menu is used to adjust the priority mprime runs at.
You should not need to change this.  You might raise the priority
if you (or your coworker) just cannot live without your screen
saver (bad idea), or if you are running some ill-behaved program
that is using CPU cycles for no good reason.

The Manual Communication menu choice should only be used if the
automatic detection of an Internet connection is not working for you.
Using this option means you have to remember to communicate with the
server every week or two (by using this same menu choice).

The Quit GIMPS menu choice is used when you no longer want this computer
to work on the GIMPS project.  You may rejoin at a later date.
If you are a PrimeNet user your unfinished work will be returned to the
server.  If you are a manual user, you need to send me email containing
your results.txt file and a note saying you are quitting.


OPTIONS MENU
------------

The CPU menu choice lets you tell the program what kind of CPU you
have and how much memory the program can use (see the earlier section
on "Setting available memory".  The program will normally figure out
the CPU type for you.  This information is used to give accurate time
estimates in the Test/Status menu choice.  It is also used, to choose
between integer based factoring (Cyrix and 486) and floating point based
factoring (Pentium and Pentium Pro).  There are also some optimizations
specific to the Pentium Pro and later CPUs in the Lucas-Lehmer primality test.
If you have an AMD K5 or Cyrix M2 chip, choose Cyrix 6x86.

The Preferences menu choice lets you control how often a line is
written to the main window and how often a line is written to
the results file.  It also lets you change how often
intermediate files (to guard against power failure and crashes)
are created.  You can control how often the program checks to
see if you are connected to the Internet.  The program polls
whenever it has new data to send to or work to get from the PrimeNet
server.  If you are low on disk space, you can select one intermediate 
file instead of two.  However, if you crash in the middle of writing
the one intermediate file, you may have to restart an exponent from
scratch.  You can also tell the program to be quiet, rather than 
beeping like crazy, if a new Mersenne prime is found.

The Torture Test choice will run a continuous self test.  This is great
for testing machines for hardware problems.  See the file stress.txt
for a more in-depth discussion of stress testing and hardware problems.

The Benchmark choice times the program on several FFT lengths.  You can
then compare your computer's speed to others list at
http://www.mersenne.org/bench.htm


POSSIBLE HARDWARE FAILURE
-------------------------

If the message "Possible hardware failure, consult the readme file."
appears in the results.txt file, then the program's error-checking has
detected a problem.  After waiting 5 minutes, the program will continue
testing from the last save file.

Could it be a software problem?  If the error is ILLEGAL SUMOUT, then 
it is very likely a true hardware error.  Although, it is impossible
to rule out an OS error corrupting part of mprime's address space.
The good news is that the program recovers very well from ILLEGAL SUMOUT
errors.

The other two errors messages, SUMINP != SUMOUT and ROUND OFF > 0.40 are
caused by one of two things:
	1)  For reasons too complicated to go into here, the program's error
	checking is not	perfect.  Some errors can be missed and some correct
	results flagged as an error.  If you get the message "Disregard last
	error..." upon continuing from the last save file, then you may have
	found the rare case where a good result was flagged as an error.
	2)  A true hardware error.

If you do not get the "Disregard last error..." message or this happens
more than once, then your machine is a good candidate for a torture test.
See the stress.txt file for more information.

Running the program on a computer with hardware problems may produce
incorrect results.  This won't hurt the GIMPS project since all results
will be double-checked.  However, you could be wasting your CPU time.
If you are getting several errors during each primality test, then I
recommend using your machine to factor Mersenne numbers or run less
strenuous distributed computing projects at
http://www.mersenne.org/projects.htm.


LUCAS-LEHMER DETAILS
--------------------

This program uses the Lucas-Lehmer primality test to see if 2**p-1 is prime.
The Lucas sequence is defined as:
	L[1] = 4
	L[n+1] = (L[n]**2 - 2) mod (2**p - 1)
2**p-1 is prime if and only if L[p-1] = 0.

This program uses a discrete weighted transform (see Mathematics of
Computation, January 1994) to square numbers in the Lucas-Lehmer sequence.


DISCLAIMER
----------

THIS PROGRAM AND INFORMATION IS PROVIDED "AS IS" WITHOUT WARRANTY OF
ANY KIND, EITHER EXPRESSED OR IMPLIED, INCLUDING BUT NOT LIMITED TO
THE IMPLIED WARRANTIES OF MERCHANTABILITY AND/OR FITNESS FOR A
PARTICULAR PURPOSE.


THANKS
------

Happy hunting and thanks for joining the search,
George Woltman
woltman@alum.mit.edu

