#define CHECK_OFTEN		0
#define MIN_K_BITS		1
#define MAX_K_BITS		49
#define MAX_C_BITS_FOR_SMALL_K	4
#define MAX_C_BITS_FOR_LARGE_K	20
#define MAX_B			2
#define MIN_N			300
#define MAX_N			12000000

/* Generate a random k of the given size */

double gen_k (int bits)
{
	double	k;
	int	i;

/* To reduce number of generic reductions due to gcd (k,c) > 1, only */
/* return small k values not divisible by 3, 5, and 7 */

	do {
		k = 1.0;
		for (i = bits; --i; )
			k = k * 2.0 + ((i == 1) ? 1 : (rand () & 1));
	} while (k > 7.0 && k < 4000000000.0 &&
		 ((unsigned long) k % 3 == 0 ||
		  (unsigned long) k % 5 == 0 ||
		  (unsigned long) k % 7 == 0));
	return (k);
}

/* Generate a random c of the given size */

int gen_c (int bits)
{
	int	c, i;

/* To reduce number of generic reductions due to gcd (k,c) > 1, only */
/* return c values not divisible by 3, 5, and 7 */

	do {
		c = 1;
		for (i = bits; --i; )
			c = c * 2 + ((i == 1) ? 1 : (rand () & 1));
	} while (c > 7 && (c % 3 == 0 || c % 5 == 0 || c % 7 == 0));
	if (rand () & 1) c = -c;
	return (c);
}

/* Generate a random n of the given size */

unsigned int gen_n ()
{
	unsigned int n, nbits;
static	unsigned int max_n_bits = 0;

	if (max_n_bits == 0)
		max_n_bits = (unsigned int) ceil (log(MAX_N-MIN_N) / log(2));
	nbits = (unsigned int) rand () % max_n_bits + 1;
	n = 1;
	while (--nbits) n = n * 2 + (rand () & 1);
	n += MIN_N;
	if (n > MAX_N) n = MAX_N;
	return (n);
}


/* Generate random data to start the test */

void gen_data (gwnum x, giant g)
{
	unsigned long i, len;

/* Generate the random number */

	len = (((unsigned long) bit_length) >> 5) + 1;
	for (i = 0; i < len; i++) {
		g->n[i] = ((unsigned long) rand() << 20) +
			  ((unsigned long) rand() << 10) +
			  (unsigned long) rand();
	}
	g->sign = len;
	specialmodg (g);
	gianttogw (g, x);
}

/* Set and print random seed */

void set_seed ()
{
	int	seed;
	char	buf[100];

	seed = time (NULL);
	srand (seed);
	sprintf (buf, "Random seed is %ld\n", seed);
	OutputBoth (buf);
}

/* Compare gwnum to a giant */

void compare (gwnum x, giant g)
{
	giant	tmp;

	tmp = popg ((PARG >> 4) + 13);
	gwtogiant (x, tmp);
	if (gcompg (g, tmp) != 0)
		OutputBoth ("Test failed.\n");
	pushg (1);
}

/* Thoroughly test the current setup */

void test_it (int square_only)
{
	gwnum	x, x2, x3, x4;
	giant	g, g2, g3, g4;
	int	i;
	double	diff, maxdiff = 0.0;
	char	buf[200];

/* Init random number generator */

//	set_seed ();

/* Alloc and init numbers */

	x = gwalloc ();
	x2 = gwalloc ();
	x3 = gwalloc ();
	x4 = gwalloc ();
	g = popg (((unsigned long) bit_length >> 4) + 20);
	g2 = popg (((unsigned long) bit_length >> 4) + 20);
	g3 = popg (((unsigned long) bit_length >> 4) + 20);
	g4 = popg (((unsigned long) bit_length >> 4) + 20);
	gen_data (x, g);
	gwcopy (x, x2); gtog (g, g2);

/* Test 50 squarings */	

	gwsetnormroutine (0, 1, 0);		/* Enable error checking */
	for (i = 0; i < 50; i++) {
		gwstartnextfft ((i & 3) == 2);	/* Test POSTFFT sometimes */
						/* Test gwsetaddin without */
						/* and with POSTFFT set */
		if ((i == 45 || i == 46) && abs (CARG) == 1)
			gwsetaddin (-31);
		gwsquare (x);
		diff = fabs (gwsuminp (x) - gwsumout (x));
		if (diff > maxdiff) maxdiff = diff;
		squareg (g);
		if ((i == 45 || i == 46) && abs (CARG) == 1) {
			iaddg (-31, g);
			gwsetaddin (0);
		}
		specialmodg (g);
		if (CHECK_OFTEN && (i & 3) != 2) compare (x, g);
	}
	if (square_only) goto done;
	if (MAXDIFF < 1e50)
		sprintf (buf, "Squares complete. MaxErr=%.8g, SumoutDiff=%.8g/%.8g(%d to 1)\n", MAXERR, maxdiff, MAXDIFF, (int) (MAXDIFF / maxdiff));
	else
		sprintf (buf, "Squares complete. MaxErr=%.10g\n", MAXERR);
	OutputBoth (buf);

/* Test mul by const */

	gwsetmulbyconst (3);
	gwsetnormroutine (0, 1, 1);
	gwsquare (x);
	gwsetnormroutine (0, 1, 0);
	diff = fabs (gwsuminp (x) - gwsumout (x));
	if (diff > maxdiff) maxdiff = diff;
	squareg (g); imulg (3, g);
	specialmodg (g);
	if (CHECK_OFTEN) compare (x, g);

/* Test square carefully */

	if (abs (CARG) == 1) gwsetaddin (-42);
	gwsquare_carefully (x);
	diff = fabs (gwsuminp (x) - gwsumout (x));
	if (diff > maxdiff) maxdiff = diff;
	squareg (g);
	if (abs (CARG) == 1) { iaddg (-42, g); gwsetaddin (0); }
	specialmodg (g);
	if (CHECK_OFTEN) compare (x, g);

/* Test gwaddquick, gwsubquick */

	gwadd3quick (x, x2, x3); gtog (g, g3); addg (g2, g3);
	gwsub3quick (x, x2, x4); gtog (g, g4); subg (g2, g4);
	if (CHECK_OFTEN) {
		specialmodg (g3); compare (x3, g3);
		specialmodg (g4); compare (x4, g4);
	}

/* Test gwadd and gwsub */

	gwadd (x, x); gwadd (x, x); gwadd (x, x); imulg (8, g);
	gwsub (x3, x); subg (g3, g);
	gwadd (x4, x); addg (g4, g);
	gwadd3 (x3, x4, x2); gtog (g3, g2); addg (g4, g2);
	gwsub3 (x3, x, x4); gtog (g3, g4); subg (g, g4);
	if (CHECK_OFTEN) {
		specialmodg (g); compare (x, g);
		specialmodg (g2); compare (x2, g2);
		specialmodg (g4); compare (x4, g4);
	}
	gwadd (x2, x); addg (g2, g);
	gwadd (x4, x); addg (g4, g);

/* Test gwaddsub */

	gwaddsub (x, x2);	// compute x+x2 and x-x2
	subg (g, g2);		// x2-x
	addg (g, g);		// x+x
	addg (g2, g);		// x+x2
	negg (g2);		// x-x2
	gwaddsub4 (x, x2, x3, x4);	// compute x+x2 and x-x2
	gtog (g, g3); addg (g2, g3);
	gtog (g, g4); subg (g2, g4);
	if (CHECK_OFTEN) {
		specialmodg (g); compare (x, g);
		specialmodg (g2); compare (x2, g2);
		specialmodg (g3); compare (x3, g3);
		specialmodg (g4); compare (x4, g4);
	}
	gwadd (x2, x); addg (g2, g);
	gwadd (x3, x); addg (g3, g);
	gwadd (x4, x); addg (g4, g);

/* Test gwaddsmall */

	gwaddsmall (x, 111);
	iaddg (111, g); specialmodg (g);
	if (CHECK_OFTEN) compare (x, g);

/* Do some multiplies to make sure that the adds and subtracts above */
/* normalized properly. */

	gwfft (x, x);
	gwfftfftmul (x, x, x);
	diff = fabs (gwsuminp (x) - gwsumout (x));
	if (diff > maxdiff) maxdiff = diff;
	squareg (g); specialmodg (g);
	if (CHECK_OFTEN) compare (x, g);

	gwfft (x, x2); gwcopy (x2, x); gwfftadd3 (x, x2, x4);
	gwfftmul (x4, x3); addg (g3, g3); mulg (g, g3); specialmodg (g3);
	diff = fabs (gwsuminp (x3) - gwsumout (x3));
	if (diff > maxdiff) maxdiff = diff;
	if (CHECK_OFTEN) compare (x3, g3);
	gwfft (x3, x4);
	gwfftfftmul (x4, x2, x);
	diff = fabs (gwsuminp (x) - gwsumout (x));
	if (diff > maxdiff) maxdiff = diff;
	mulg (g3, g); specialmodg (g);
	if (CHECK_OFTEN) compare (x, g);

/* Do the final compare */

done:	if (!CHECK_OFTEN) compare (x, g);
	if (GWERROR) OutputBoth ("GWERROR set during calculations.\n");

/* Print final stats */

	if (maxdiff > MAXDIFF) OutputBoth ("Sumout failed during test.\n");
	if (MAXDIFF < 1e50)
		sprintf (buf, "Test complete. MaxErr=%.8g, SumoutDiff=%.8g/%.8g(%d to 1)\n", MAXERR, maxdiff, MAXDIFF, (int) (MAXDIFF / maxdiff));
	else
		sprintf (buf, "Test complete. MaxErr=%.10g\n", MAXERR);
	OutputBoth (buf);
	OutputBoth ("\n");
}


void test_methodically (void)
{
	int	kbits, cbits;
	double	k;
	unsigned int b;
	int	c;
	unsigned int n;

// Test general mods
// Test rational mods
// Test near fftlen max and test near fft min

	b = 2;
//	for k = 1 bit to 48 bits
//		for c = 1 bit to +/-24 bits
//			for n = ? to ? to test each fftlen

	set_seed ();
	k = gen_k (kbits);
	c = gen_c (cbits);
	gwsetup (k, b, n, c, 0);
// print number as string
// print fftlen
	test_it (1);
	gwdone ();
}

void test_randomly (void)
{
	int	kbits, cbits;
	double	k;
	unsigned int b;
	int	c;
	unsigned int n;
	char	buf[100], fft_desc[100];

	set_seed ();
	b = 2;
	for ( ; ; ) {
		if (MAX_K_BITS == MIN_K_BITS) kbits = 0;
		else kbits = (unsigned int) rand () % (MAX_K_BITS - MIN_K_BITS + 1);
		kbits += MIN_K_BITS;
		if (kbits <= 20)
			cbits = (unsigned int) rand () % MAX_C_BITS_FOR_SMALL_K + 1;
		else
			cbits = (unsigned int) rand () % MAX_C_BITS_FOR_LARGE_K + 1;
		k = gen_k (kbits);
		c = gen_c (cbits);
		n = gen_n ();
		switch (rand () % 5) {
		case 0:	CPU_L2_CACHE_SIZE = 128; break;
		case 1:	CPU_L2_CACHE_SIZE = 256; break;
		case 2:	CPU_L2_CACHE_SIZE = 512; break;
		case 3:	CPU_L2_CACHE_SIZE = 1024; break;
		case 4:	CPU_L2_CACHE_SIZE = -1; break;
		}
		gwsetup (k, b, n, c, 0);
		sprintf (buf, "Starting QA run on %s, kbits=%ld, cbits=%ld\n",
			 gwmodulo_as_string (), kbits, cbits);
		OutputBoth (buf);
		gwfft_description (fft_desc);
		if (FFTLEN < 600000)
			sprintf (buf, "Using %s, bits_per_word=%.5g/%.5g\n",
				 fft_desc, virtual_bits_per_word (),
				 fft_max_bits_per_word);
		else
			sprintf (buf, "Using %s, bits_per_word=%.5g/%.5g, L2_cache_size=%d\n",
				 fft_desc, virtual_bits_per_word (),
				 fft_max_bits_per_word, CPU_L2_CACHE_SIZE);
		OutputBoth (buf);
		test_it (0);
		gwdone ();
	}
}
