
Welcome to the Great Internet Mersenne Prime Search!

To use this program you must agree to the prize rules at
http://www.mersenne.org/prize.htm

In case you ever forget, the URL is http://www.mersenne.org/prime.htm.
My email address is woltman@alum.mit.edu.  For networking questions,
contact Scott Kurowski's technical support team at primenet@entropia.com.


FILE LIST
---------

readme.txt	This file.
ntsetup.exe	This program is used to configure ntprime.exe options.
ntprime.exe	This is the NT Service that runs the Lucas-Lehmer tests
		to find Mersenne primes.
whatsnew.txt	A list of new features in ntprime.exe.
prime.ini	A file containing your preferences.  The menu choices
		and dialog boxes are used to change your preferences.
local.ini	Like prime.ini, this file contains more preferences.
		The reason there are two files is discussed later.
worktodo.ini	A list of exponents the program will be factoring
		and/or Lucas-Lehmer testing.
results.txt	Ntprime.exe writes its results to this file.
prime.log	A text file listing all messages that have been sent
		to the PrimeNet server.
prime.spl	A binary file of messages that have not yet been sent to
		the PrimeNet server.
pNNNNNNN &	Intermediate files produced by ntprime.exe to resume
qNNNNNNN	computation where it left off.
mNNNNNNN	Intermediate files produced during P-1 factoring.


WHAT IS THIS PROGRAM?
---------------------

This program is used to find Mersenne Prime numbers.  See
http://www.utm.edu/research/primes/mersenne.shtml for a good
description of Mersenne primes.  Mersenne numbers can be proved
composite (not prime) by either finding a factor or by running
a Lucas-Lehmer primality test.


INSTRUCTIONS
------------

If you are "upgrading" from the Prime95.exe program, you can skip
to the section on upgrading from prime95.exe.

There are two ways to use this program.  The automatic way uses
a central server, which we call the PrimeNet server, to get work to do
and report your results.  Anyone with Internet access, including AOL,
CompuServe, and dial-up ISP users should use this method.  You do not
need a permanent connection to the Internet.

The second method is the manual method.  It requires a little more work
and monitoring.  I recommend this for computers with no Internet access
or with some kind of firewall problem that prevents the automatic method
from working.

If you are running this program at your place of employment, you should
first GET PERMISSION from your network administrator or boss.  This is
especially true if you are installing the software on several machines.
Some companies are reluctant to run any software they are not familiar with.

Note:  The instructions below are for Windows NT 3.51, Windows NT 4.0 may
differ slightly.


INSTRUCTIONS FOR UPGRADING FROM AN EARLIER NT SERVICE VERSION
-------------------------------------------------------------

1)  Stop the old service.
2)  Type "NTPrime -remove" to remove the old service.
3)  Unzip the new version on top of the old version.
4)  Type "NTPrime -install" to install the new version.
5)  Start the new service.


INSTRUCTIONS FOR THE AUTOMATIC METHOD
-------------------------------------

1)  Download and unzip winnt.zip.  You've probably done this already
    since you are reading this file.
2)  Connect to the Internet.
3)  Run ntsetup.exe.  You will see 3 dialog boxes:
3a) In the first dialog box, enter your name and email address.
    Optionally enter a user ID, password, and computer ID.  If you are using
    several computers, use the same user ID and password but a unique computer
    ID on each machine.  An easy-to-remember user ID will be helpful if you
    plan to visit the PrimeNet server's web page to view reports on your
    progress.  If you do not enter a user ID or if you pick a user ID that
    is already in use, then the server will assign a user ID.
3b) In the second dialog box, fill in roughly how many hours a day you leave
    your computer running.  Click OK.
3c) In the third dialog box, leave the "Use Primenet..." checkbox
    checked.  Do not turn this checkbox off even if you disconnect
    from the Internet.   Check the "Use a dial-up..." checkbox
    if you use a modem to connect to the Internet.  Note that ntprime
    will not dial-up to connect to the Internet, rather it waits
    for a time when you are already connected to contact the server.
    Click OK.  Ntprime will now contact the PrimeNet server to get some
    work for your computer to do.
4)  Try contacting the PrimeNet server by choosing
    Preferences/Manual Communication.  Hit OK.
4a) If a proxy server is causing connection troubles, see the
    later section on "SETTING UP A PROXY SERVER".
4b) If you cannot contact the PrimeNet server,
    you will have to use the manual method described below.
5)  Exit ntsetup.exe.
6)  Run "ntprime -install" to install the Service.
7)  From Services in the Control Panel, select "Prime Service" and click
    Start.
8)  Clicking Stop will write intermediate results to a file, and then stop the
    service.  Clicking Start will restart the service where it left off.


MANUAL METHOD INSTRUCTIONS
--------------------------

1)  Use the Web (http://entropia.com/ips/manualtests.html) to create
    a userid for yourself and to get a set of exponents to work on.
    Alternatively, you can get some work from
    http://www.mersenne.org/range2.htm.  Copy these exponents to a
    file called worktodo.ini
2)  Run ntsetup.exe.  You will see 3 dialog boxes:
2a) In the first dialog box, enter your name and email address.  Click OK.
2b) In the second dialog box, fill in roughly how many hours a day you leave
    your computer running.  Click OK.
2c) In the third dialog box, uncheck "Use PrimeNet to get work and report
    results", click OK.
3)  Exit ntsetup.exe.
4)  Run "ntprime -install" to install the Service.
5)  From Services in the Control Panel, select "Prime Service" and click
    Start.  Clicking Stop will write intermediate results to a file, and
    then stop the service.  Clicking Start will restart the service where
    it left off.
6)  Once a month or when done with your exponents, use the web pages
    again to send the file "results.txt" to the PrimeNet server.
    It is important to do this so the exponents you are testing are
    not reassigned to someone else.


INSTRUCTIONS FOR UPGRADING FROM VERSION 15 OF PRIME95.EXE
---------------------------------------------------------

1)  Download and unzip winnt.zip in the same directory as prime95.exe.
2)  Run "ntprime -install" to install the Service.
3)  From Services in the Control Panel, select "Prime Service" and click
    Start.
4)  Clicking Stop will write intermediate results to a file, and then stop the
    service.  Clicking Start will restart the service where it left off.


NOTES ON NT SECURITY
--------------------

Only members of the Administrators group can install new services (I think
installing new services requires the SE_LOAD_DRIVER_NAME privilage).  As it
is currently set up, only Administrators can start, stop, or configure the
Prime service, although this may change in the future (Power Users can start
and stop, maybe?).

If you're interested in keeping a secure system, you probably already know
the following (and if you're not interested, you can ignore the following):

The ntprime.exe file should give Full Control to Adminstrators and SYSTEM,
and Read access to Everyone.  It should be owned by Adminstrators.

It is probably a good idea to to create a new account for the Prime Service
to run in (as apposed to letting the service run in the omnipotent
LocalSystem account).  In order to do this, create a new account called
GuestService (or whatever) and give it a password you can remember for the
next few minutes.  GuestService should be a member of the group Guests.
Next, give GuestService the "Log on as a service" privilege.  (In the User
Manager:  select Policies | User Rights...; Show Advanced User Rights;
next to "Right:", select "Log on as a service"; Add...; type GuestService
under "Add Names:"; OK; OK.)  Next, from Services in the Control Panel,
select "Prime Service" and then click "Startup...".  Click "This Account"
and then type GuestService, and then enter the password you chose.  If the
service was stopped, try starting it to make sure you typed the password
correctly.


OTHER NOTES
-----------

Let ntprime.exe run at all times.  It runs at the lowest possible priority,
making use of all your idle CPU cycles.  It should not interfere with your
normal work.  Let the program run overnight and on weekends.  Never turn
your computer off.  Turn off your monitor to conserve energy.  NOTE:
Running your computer non-stop could increase your electric bill by $30
per year or more.

In order for new changes to take effect after running ntsetup, you must
restart (Stop and then Start) the service if it is already running.

To start 2 services on a dual-processor machine, create a second
directory and follow all the instructions above, except that before
running "ntprime -install" you first add the lines
"ServiceName=NTPrimeService2" and
"DisplayName=Prime Service 2" to the prime.ini file.

You can have ntprime contact the PrimeNet server in a batch file.  Execute
"ntprime -c" to have ntprime contact the server and then exit.

You can get a report of your PrimeNet server account status on the web
(http://entropia.com/primenet/status.shtml).  Type your UserID and
password into the web form, and click "Get Account Report".  You can
get your UserID and password from the Test/User Information
dialog box.

Information about running a local PrimeNet server is also available.
See http://entropia.com/primenet.  If you have any questions about the
PrimeNet server, you can send e-mail to primenet@entropia.com.

A line is written to the results file every time the service is stopped
so that you can monitor the program's progress.

The first time you run this program it will perform an hour long
self-test to make sure the Lucas-Lehmer code is running properly
on your system.

If you have overclocked your machine, I highly recommend running the
torture test for a couple of days.  The longer you run the torture test
the greater the chance that you will uncover an error caused by
overheating.

Depending on the exponent being tested, the program may decide that it
would be wise to invest some time checking for small factors before
running a Lucas-Lehmer test.  Furthermore, the program may start factoring
exponents before a previous Lucas-Lehmer test completes.  This is normal!
The program will resume the Lucas-Lehmer test when the factoring
completes.

You can compare your computer's speed with other users by checking the
site http://www.mersenne.org/bench.htm.  

It can take many CPU days to test a large Mersenne number.  The service
can be safely interrupted to write intermediate results to disk.  This
program also saves intermediate results to disk every 30 minutes
in case there is a power failure.

You can remove the service by running ntprime -remove.

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

Create a file in the same local folder as Ntprime.exe, called "primenet.ini".
Add these text lines, substituting the appropriate proxy server URL and port:

[PrimeNet Proxy]
ProxyHost=http://proxy.megacorp.com:8080

If the proxy is secured by a userid/password, add the following two text
lines, substituting the appropriate values:

ProxyUser=whatever_the_userid_is
ProxyPass=whatever_the_password_is

After the first time PrimeNet is contacted through a secure proxy, the proxy 
password is encoded and a new parameter 'ProxyMask=1' set.  To change the 
password, simply change the ProxyPass= value, and either delete the ProxyMask 
text line or set ProxyMask=0.


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


RUNNING NTPRIME ON SEVERAL COMPUTERS
------------------------------------

The easiest way to do this is to first set up ntprime on one computer.
Next copy all the files to the second computer.  Delete the local.ini
file and worktodo.ini files.  These files contain information that
is specific to the first computer.  Start ntprime on the second
computer and optionally use Test/User Information to
give the second computer a unique computer ID.  Repeat this process
for all the computers you wish to run ntprime on.

If you do not follow the instruction above, be sure you use 
Test/User Information to give each computer the same
userid and password.  Failure to do this will result in all 
your work being "credited" to different user IDs.


PROGRAM OUTPUT
--------------

The results file will include lines that look like:

Factored M400037 through 17517*2^32 (pass 3 of 16).
	This means ntprime is in the third pass of a 16 pass process to
	find a small factor of 2^400037-1.  
Iteration: 941400 / 1667747.
	This means ntprime just finished the 941400th iteration of a
	Lucas-Lehmer primality test.  The program must execute 1667747
	iterations to complete the primality test.
M2645701 has a factor: 13412891051374103
	This means to 2^2645701-1 is not prime.  It is divisible
	by 13412891051374103.
M2123027 no factor to 2^57, WS5: 14780E25
	This means 2^2123027-1 has no factors less than 2^57.  The Mersenne
	number may or may not be prime.  A Lucas-Lehmer test is needed
	to determine the primality of the Mersenne number.  WS5 is
	the program version number.  14780E25 is a checksum to guard
	against transmission errors.
M1992031 is not prime. Res64: 6549369F4962ADE0. WS5: B253EF24
	This means 2^1992031-1 is not prime - a Lucas-Lehmer test says so.
	The last 64 bits of the last number in the Lucas-Lehmer sequence
	is 6549369F4962ADE0.  At some future date, another person using
	another program will run the same Lucas-Lehmer test.  If his last 64
	bits match yours, then both machines ran flawlessly and it is 100%
	certain that the Mersenne number is composite.  WS5 is the program
	version number.  B253EF24 is a checksum to guard against
	transmission errors.
M11213 is prime! WS5: 579A579A
	This means 2^11213-1 is a Mersenne prime!  WS5 is the program
	version number.  579A579A is a checksum to guard against
	transmission errors.


NTSETUP TEST MENU
-----------------

The PrimeNet menu choice lets you configure the type of work you get
from the server.  The "Use PrimeNet..." option can be turned on to
switch from the manual method to the automatic method.
The "Request whatever..." box should be left checked.  However, if you
are running a slow computer and don't mind waiting several months for
a single Lucas-Lehmer test to complete OR you are running a faster
computer and would rather do factoring, then uncheck this box and
choose a different type of work to do.  The "Always have at least
this many days of work queued up" should be changed based on how often
you connect to the Internet.  As long as you connect at least once in the
given time period, ntprime will have an uninterrupted stream of work.
However, the program will not checkout more than 20 exponents no matter
what this value is set at.

The P-1 choice lets you factor Mersenne numbers using the P-1 method of
factoring.  There is presently no web site which tells you how much
P-1 factoring has already been done on exponents.  You can also edit
the worktodo.ini file directly.  For example:
	Pminus1=751001,1000000,0,0,0
The first value is the exponent.  The second value is bound #1.  The
third value is bound #2.  The fourth value is 0 for 2^N-1 factoring,
1 for 2^N+1 factoring.  The fifth value is no longer used.

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

The Quit GIMPS menu choice is used when you no longer want this computer
to work on the GIMPS project.  You may rejoin at a later date.
If you are a PrimeNet user your unfinished work will be returned to the
server.  If you are a manual user, you need to send me email containing
your results.txt file and a note saying you are quitting.

The User Information menu choice lets you change your name and email
address.  Your name will be used in credits and "top producers"
web pages.  Your email address will be used to send you warnings if
an exponent you have reserved is about to expire and to send you a
newsletter once every few months.  Select the no email checkbox to
stop the server from sending expiration warnings.  The userid you
use will appear in the status pages on the PrimeNet server.  You can
also set the ComputerID field if you have several computers and want
to keep track of which computers produced which results.

The Vacation/Holiday menu choice lets ntprime update the expected completion
dates on the PrimeNet server.  This will prevent one of your exponents
from being reassigned while you're gone.  Also, if you are taking a
long vacation, ntprime will get extra exponents to test.  For example,
if you are leaving immediately for a 6 week vacation, connect to the
Internet and choose 42 days.  If you are leaving in two days for a 6
week vacation, then set choose 44 days.  Primenet will send the new
completion dates whenever you connect to the Internet over the next 2 days.

The Status menu choice will tell you what exponents you are working on.
It will also estimate how long that will take and
your chances of finding a new Mersenne prime.


NTSETUP OPTIONS MENU
--------------------

The CPU menu choice lets you tell the program what kind of CPU you
have and how much memory the program can use (see the earlier section
on "Setting available memory".  The program will normally figure out
the CPU type for you.  This information is used to give accurate time
estimates in the Test/Status menu choice.  It is also used, to choose
between integer based factoring (Cyrix and 486) and floating point based
factoring (Pentium and Pentium Pro).  There are also some optimizations
specific to the Pentium Pro and later CPUs in the Lucas-Lehmer primality test.
If you have an AMD K5 or Cyrix M2 chip, choose Cyrix 6x86.

The Preferences menu choice lets you control how often a line is written
to the results file.  It also lets you change how often
intermediate files (to guard against power failure an crashes)
are created.  You can control how often the program checks to
see if you are connected to the Internet.  The program polls
whenever it has new data to send to or work to get from the PrimeNet
server.  If you are low on disk space, you can select one intermediate 
file instead of two.  However, if you crash in the middle of writing
the one intermediate file, you may have to restart an exponent from
scratch.  You can also tell the program to be quiet, rather than 
beeping like crazy, if a new Mersenne prime is found.

The Manual Communication menu choice allows you to control when you
contact the PrimeNet server.  Check the first checkbox if you never
want ntprime.exe to automaticly contact the server.  This is not
recommended as using this option means you have to remember to
communicate with the server every week or two.  The second checkbox
tells ntsetup.exe to contact the server right now.

The Priority menu is used to adjust the priority prime95 runs at.
You should not need to change this.  You might raise the priority
if you (or your coworker) just cannot live without your screen
saver (bad idea), or if you are running some ill-behaved program
that is using CPU cycles for no good reason.


POSSIBLE HARDWARE FAILURE
-------------------------

If the message, "Possible hardware failure, consult the readme file.",
appears in the results file, then ntprime's error-checking has
detected a problem.  Ntprime will continue from the last save file.
If you do not get the message, "Disregard last error...", then the
problem is not reproducible - a definite sign of hardware problems.

How can this be when none of your other programs have problems?  The answer
is that ntprime stresses your machine more than any other program you
run.  The operating system usually shuts down the floating-point unit
when no programs are using it.  Ntprime continuously uses the FPU, consuming
more electricity and generating more heat.  If the CPU is not properly cooled,
errors can occur.  Ntprime also constantly accesses main memory - up to
60MB per second.  This constant activity will detect memory problems that
other programs do not.  This is why Cray Research has used a program similar
to this one as part of its supercomputer diagnostics package for over a decade.

Could it be a software problem?  If the error is ILLEGAL SUMOUT, then
there is a good chance that this is a software problem.  A device driver
or VxD may not be saving and restoring CPU state correctly.  The good news
is that ntprime recovers very well from ILLEGAL SUMOUT errors.
Try seeing if the problem occurs only when a specific device is active
or a specific program is running.

How can you track down the hardware problem?  Unfortunately, this is not
easy.  To see if your CPU is overheating, run ntprime for several hours.
Open the box.  Is the CPU too hot to touch?  If so, a heat sink or
CPU fan should solve the problem.  Memory problems are not as easy to
diagnose.  My only advice is to try swapping memory SIMMs with a coworker's
or friend's machine.  If the errors go away, then you can be confidant
that the original problems were memory related.

What can you do if you are unwilling or unable to find the hardware problem?
If you are only getting an error once in a while, then your results are
probably OK.  The error-checking code is not infallible, so your results
will need to be double-checked.  If you are getting several errors during
each primality test, then I would recommend using your machine to factor
Mersenne numbers.


DETAILS
-------

This program uses the Lucas-Lehmer method to test if 2**p-1 is prime.  The
Lucas sequence is defined as:
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
