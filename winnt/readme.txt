Welcome to the Great Internet Mersenne Prime Search!

Remember:  If you are running this program at your place of employment, you
must first GET PERMISSION from your network administrator or boss.
Some companies are reluctant to run any software they are not familiar with.
Failing to get permission is often grounds for dismissal.

This NT Service version is only for advanced Windows NT/2000/XP users.
The normal client, prime95.exe, can run as a GUI NT service simply by
selecting the Options/Start at Bootup option.   This service version
is provided for those administrators wary of allowing a GUI based service
run in the Local System account.

INSTRUCTIONS
------------

Install and configure prime95:

1) You must also install and run the normal GUI-based GIMPS client,
   prime95.exe, to configure this service.  Install this prior to
   installing the service.
2) Follow the prime95.exe instructions to configure the client.
3) Make sure the Options/Start at Bootup option is not set.
4) Exit prime95.exe

Now install the service:

1) Unzip the ntprime.exe program in the same directory as prime95.exe.
2) Run "ntprime -install" to install the Service.
3) From Services in the Control Panel, select "Prime Service" and click Start.

To configure the service after it is installed:

1) Stop the service.
2) Run Prime95, make your configuration changes.
3) Exit prime95.
4) Start the service.


NOTES ON NT SECURITY
--------------------

These notes are old and provided by a GIMPS member.  I am not knowlegable
enough in NT administration to comment on their accuracy.

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

To start 2 services on a dual-processor machine, create a second
directory and follow all the instructions above, except that before
running "ntprime -install" you first add the lines
"ServiceName=NTPrimeService2" and
"DisplayName=Prime Service 2" to the prime.ini file.

You can have ntprime contact the PrimeNet server in a batch file.  Execute
"ntprime -c" to have ntprime contact the server and then exit.

A line is written to the results file every time the service is stopped
so that you can monitor the program's progress.

You can remove the service by running ntprime -remove.


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
