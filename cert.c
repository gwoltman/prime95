/*----------------------------------------------------------------------
| Copyright 2020-2026 Mersenne Research, Inc.  All rights reserved
|
| This file contains routines to certify a PRP proof.
+---------------------------------------------------------------------*/

/* Cert state to be written to / read from certify save files */

struct cert_state {
	int	thread_num;		/* Copy of thread_num passed to cert.  Allows subroutines to output error messages. */
	gwhandle gwdata;
	gwnum	x;			/* The current value in our squaring loop */
	unsigned long counter;		/* Current "iteration" counter */
	unsigned long units_bit;	/* For shifting FFT data -- allows more robust double-checking */
	unsigned long error_count;	/* Count of any errors that have occurred during certify test */
	int	state;			/* Gerbicz or double-check state, uses the PRP_STATE_* definitions */
	gwnum	alt_x;			/* When double-checking, this is the alternate x value */
					/* When Gerbicz checking, this is the comparison checksum value being calculated */
	unsigned long alt_units_bit;	/* When double-checking, this is the alternate shift count */
	gwnum	u0;			/* Saved first value of the Gerbicz checksum function */
	gwnum	d;			/* Last computed value of the Gerbicz checksum function */
	unsigned long L;		/* Iterations between multiplies in computing Gerbicz checksum */
	unsigned long start_counter;	/* Counter at start of current Gerbicz or double-check block */
	unsigned long next_mul_counter;	/* Counter when next Gerbicz multiply takes place */
	unsigned long end_counter;	/* Counter when current Gerbicz or double-check block ends */
};

/* Write intermediate certify results to a file */
/* The cert save file format is: */
/*	u32		magic number  (different for ll, p-1, prp, tf, ecm) */
/*	u32		version number */
/*	double		pct complete */
/*	char(11)	stage */
/*	char(1)		pad */
/*	u32		checksum of following data */
/*	u32		error_count */
/*	u32		iteration counter */
/*	u32		shift_count */
/*	u32		state (version number >= 2) */
/*	u32		alternate shift count (version number >= 2) */
/*	u32		L - iterations between Gerbicz multiplies (version number >= 2) */
/*	u32		error-checking start counter (version number >= 2) */
/*	u32		error-checking next Gerbicz multiply counter (version number >= 2) */
/*	u32		error-checking end counter (version number >= 2) */
/*	gwnum		FFT data for x (u32 len, array u32s) */
/*	gwnum		FFT data for alt_x (u32 len, array u32s) (version number >= 2) */
/*	gwnum		FFT data for u0 (u32 len, array u32s) (version number >= 2) */
/*	gwnum		FFT data for d (u32 len, array u32s) (version number >= 2) */

#define CERT_MAGICNUM		0x8f729ab1
#define CERT_VERSION		2

int writeCertSaveFile (			/* Returns TRUE if successful */
	struct cert_state *cs,		/* Certification state structure to read and fill in */
	writeSaveFileState *write_save_file_state,
	struct work_unit *w)		/* Work unit */
{
	int	fd;
	uint32_t sum = 0;
	char	buf[512];

/* Now save to the intermediate file */

	fd = openWriteSaveFile (write_save_file_state);
	if (fd < 0) {
		sprintf (buf, "Unable to create cert save file: %s\n", write_save_file_state->base_filename);
		OutputBoth (cs->thread_num, buf);
		OutputBothErrno (cs->thread_num);
		return (FALSE);
	}

	if (!write_header (fd, CERT_MAGICNUM, CERT_VERSION, w)) goto writeerr;
	if (!write_long (fd, cs->error_count, &sum)) goto writeerr;
	if (!write_long (fd, cs->counter, &sum)) goto writeerr;
	if (!write_long (fd, cs->units_bit, &sum)) goto writeerr;
	if (!write_long (fd, cs->state, &sum)) goto writeerr;
	if (!write_long (fd, cs->alt_units_bit, &sum)) goto writeerr;
	if (!write_long (fd, cs->L, &sum)) goto writeerr;
	if (!write_long (fd, cs->start_counter, &sum)) goto writeerr;
	if (!write_long (fd, cs->next_mul_counter, &sum)) goto writeerr;
	if (!write_long (fd, cs->end_counter, &sum)) goto writeerr;
	if (!write_gwnum (fd, &cs->gwdata, cs->x, &sum)) {
		sprintf (buf, "Error writing FFT data named x to cert save file %s\n", write_save_file_state->base_filename);
		OutputBoth (cs->thread_num, buf);
		OutputBothErrno (cs->thread_num);
		goto err;
	}

	if (cs->state != PRP_STATE_GERB_MID_BLOCK && cs->state != PRP_STATE_GERB_MID_BLOCK_MULT) {
		if (!write_gwnum (fd, &cs->gwdata, cs->alt_x, &sum)) {
			sprintf (buf, "Error writing FFT data named alt_x, state %d, to cert save file %s\n",
				 cs->state, write_save_file_state->base_filename);
			OutputBoth (cs->thread_num, buf);
			OutputBothErrno (cs->thread_num);
			goto err;
		}
	}

	if (cs->state != PRP_STATE_DCHK_PASS1 && cs->state != PRP_STATE_DCHK_PASS2 &&
	    cs->state != PRP_STATE_GERB_START_BLOCK && cs->state != PRP_STATE_GERB_FINAL_MULT) {
		if (!write_gwnum (fd, &cs->gwdata, cs->u0, &sum)) {
			sprintf (buf, "Error writing FFT data named u0, state %d, to cert save file %s\n",
				 cs->state, write_save_file_state->base_filename);
			OutputBoth (cs->thread_num, buf);
			OutputBothErrno (cs->thread_num);
			goto err;
		}
	}

	if (cs->state != PRP_STATE_DCHK_PASS1 && cs->state != PRP_STATE_DCHK_PASS2 && cs->state != PRP_STATE_GERB_START_BLOCK) {
		if (!write_gwnum (fd, &cs->gwdata, cs->d, &sum)) {
			sprintf (buf, "Error writing FFT data named d, state %d, to cert save file %s\n",
				 cs->state, write_save_file_state->base_filename);
			OutputBoth (cs->thread_num, buf);
			OutputBothErrno (cs->thread_num);
			goto err;
		}
	}

	if (!write_footer (fd, w)) goto writeerr;
	if (!write_checksum (fd, sum)) goto writeerr;

	closeWriteSaveFile (write_save_file_state, fd);
	return (TRUE);

/* An error occurred.  Delete the current file. */

writeerr:
	sprintf (buf, WRITEFILEERR, write_save_file_state->base_filename);
	OutputBoth (cs->thread_num, buf);
	OutputBothErrno (cs->thread_num);
err:	deleteWriteSaveFile (write_save_file_state, fd);
	return (FALSE);
}

/* Read the data portion of an intermediate cert save file */

int readCertSaveFile (			/* Returns TRUE if succsessful */
	struct cert_state *cs,		/* Certification state structure to read and fill in */
	char	*filename,		/* Save file name */
	struct work_unit *w)		/* Work unit */
{
	int	fd;
	uint32_t version, sum = 0, filesum;

	// Open the save file
	fd = _open (filename, _O_BINARY | _O_RDONLY);
	if (fd <= 0) return (FALSE);

	// Read the header
	if (!read_magicnum (fd, CERT_MAGICNUM)) goto err;
	if (!read_header (fd, &version, w, &filesum)) goto err;
	if (version == 0 || version > CERT_VERSION) goto err;

	// Read fields that are in all versions of the save file
	if (!read_long (fd, &cs->error_count, &sum)) goto err;
	if (!read_long (fd, &cs->counter, &sum)) goto err;
	if (!read_long (fd, &cs->units_bit, &sum)) goto err;

	// Version 1 save files did no Gerbicz error checking.  Start a Gerbicz block from the save file's x value.
	if (version == 1) {
		cs->state = PRP_STATE_GERB_START_BLOCK;
		if (!read_gwnum (fd, &cs->gwdata, cs->x, &sum)) goto err;
		gwcopy (&cs->gwdata, cs->x, cs->alt_x);
		cs->alt_units_bit = cs->units_bit;
	}

	// Read in fields introduced in version 2 save files
	else {
		unsigned long state;
		if (!read_long (fd, &state, &sum)) goto err;
		cs->state = (int) state;
		if (cs->state != PRP_STATE_DCHK_PASS1 && cs->state != PRP_STATE_DCHK_PASS2 &&
		    (cs->state < PRP_STATE_GERB_START_BLOCK || cs->state > PRP_STATE_GERB_FINAL_MULT)) goto err;
		if (!read_long (fd, &cs->alt_units_bit, &sum)) goto err;
		if (!read_long (fd, &cs->L, &sum)) goto err;
		if (!read_long (fd, &cs->start_counter, &sum)) goto err;
		if (!read_long (fd, &cs->next_mul_counter, &sum)) goto err;
		if (!read_long (fd, &cs->end_counter, &sum)) goto err;

		// All cert states wrote an x value
		if (!read_gwnum (fd, &cs->gwdata, cs->x, &sum)) goto err;

		// Most cert states wrote an alt_x value
		if (cs->state != PRP_STATE_GERB_MID_BLOCK && cs->state != PRP_STATE_GERB_MID_BLOCK_MULT) {
			if (!read_gwnum (fd, &cs->gwdata, cs->alt_x, &sum)) goto err;
		}

		// Most Gerbicz states wrote a u0 value
		if (cs->state != PRP_STATE_DCHK_PASS1 && cs->state != PRP_STATE_DCHK_PASS2 &&
		    cs->state != PRP_STATE_GERB_START_BLOCK && cs->state != PRP_STATE_GERB_FINAL_MULT) {
			if (!read_gwnum (fd, &cs->gwdata, cs->u0, &sum)) goto err;
		}

		// Most Gerbicz states wrote a d value
		if (cs->state != PRP_STATE_DCHK_PASS1 && cs->state != PRP_STATE_DCHK_PASS2 && cs->state != PRP_STATE_GERB_START_BLOCK) {
			if (!read_gwnum (fd, &cs->gwdata, cs->d, &sum)) goto err;
		}
	}

	// Validate checksum and return
	if (filesum != sum) goto err;
	if (!read_footer (fd, w)) goto err;
	_close (fd);
	return (TRUE);
err:	_close (fd);
	return (FALSE);
}

/* Do a PRP proof certification */

int cert (
	int	thread_num,		/* Worker number */
	struct PriorityInfo *sp_info,	/* SetPriority information */
	struct work_unit *w,		/* Worktodo entry */
	int	pass)			/* PrimeContinue pass */
{
	struct cert_state cs;
	int	first_iter_msg, res, stop_reason;
	int	echk, near_fft_limit, sleep5;
	unsigned long iters;
	double	timers[2];
	double	inverse_explen;
	double	reallyminerr = 1.0;
	double	reallymaxerr = 0.0;
	unsigned long last_counter = 0xFFFFFFFF;/* Iteration of last error */
	int	maxerr_recovery_mode = 0;	/* Big roundoff err rerun */
	double	last_maxerr = 0.0;
	long	restart_counter = -1;		/* On a restart, this specifies how far back to rollback save files */
	unsigned long restart_error_count = 0;	/* On a restart, use this error count rather than the one from a save file */
	readSaveFileState read_save_file_state; /* Manage savefile names during reading */
	writeSaveFileState write_save_file_state; /* Manage savefile names during writing */
	char	filename[32];
	char	buf[1000], JSONbuf[4000], fft_desc[200];
	double	allowable_maxerr, output_adjustment, title_adjustment;
	char	string_rep[80];
	int	error_count_messages;
	hash256_t hash;

/* We only expect certifications on Mersenne numbers from the primenet server */

	if (w->k != 1.0 || w->b != 2 || w->c != -1) {
		// If primenet ever supports non-Mersennes, we need to disable the units_bit shift and probably
		// use a safety margin of 0.25 as we do not have any error recovery from a roundoff error
		OutputStr (thread_num, "Certification of non-Mersenne numbers not supported\n");
		memset (&cs, 0, sizeof (cs));
		gwinit (&cs.gwdata);
		goto abandon_work;
	}

/* Init the write save file state.  This remembers which save files are Gerbicz-checked.  Do this initialization */
/* before the restart for errors so that error recovery does not destroy the write save file state. */

	tempFileName (w, filename);
	writeSaveFileStateInit (&write_save_file_state, filename, NUM_JACOBI_BACKUP_FILES);

/* Init cert state */

begin:	memset (&cs, 0, sizeof (cs));
	cs.thread_num = thread_num;
	gwinit (&cs.gwdata);

/* Init the FFT code for squaring modulo k*b^n+c */

	if (IniGetInt (INI_FILE, "UseLargePages", 0)) gwset_use_large_pages (&cs.gwdata);
	if (HYPERTHREAD_LL) sp_info->normal_work_hyperthreading = TRUE, gwset_will_hyperthread (&cs.gwdata, 2);
	gwset_bench_cores (&cs.gwdata, HW_NUM_CORES);
	gwset_bench_workers (&cs.gwdata, NUM_WORKERS);
	if (ERRCHK) gwset_will_error_check (&cs.gwdata);
	else gwset_will_error_check_near_limit (&cs.gwdata);
	gwset_num_threads (&cs.gwdata, get_worker_num_threads (thread_num, HYPERTHREAD_LL));
	gwset_thread_callback (&cs.gwdata, SetAuxThreadPriority);
	gwset_thread_callback_data (&cs.gwdata, sp_info);
	gwset_minimum_fftlen (&cs.gwdata, w->minimum_fftlen);
	gwset_safety_margin (&cs.gwdata, IniGetFloat (INI_FILE, "CertificationSafetyMargin", 0.0));
	gwset_use_spin_wait (&cs.gwdata, IniGetInt (INI_FILE, "SpinWait", 0));
	res = gwsetup (&cs.gwdata, w->k, w->b, w->n, w->c);

/* If we were unable to init the FFT code, then print an error message and return an error code. */

	if (res) {
		char	string_rep[80];
		gw_as_string (string_rep, w->k, w->b, w->n, w->c);
		sprintf (buf, "Certification cannot initialize FFT code for %s, errcode=%d\n", string_rep, res);
		OutputBoth (thread_num, buf);
		gwerror_text (&cs.gwdata, res, buf, sizeof (buf) - 1);
		strcat (buf, "\n");
		OutputBoth (thread_num, buf);
		goto abandon_work;
	}

/* Record the amount of memory being used by this thread.  Assume Gerbicz error-checking uses double even */
/* though it really doesn't change the working set all that much. */

	set_memory_usage (thread_num, 0, cvt_gwnums_to_mem (&cs.gwdata, 1) * 2);

/* Allocate memory for the certification and Gerbicz error checking */

	cs.x = gwalloc (&cs.gwdata);
	if (cs.x == NULL) {
		OutputStr (thread_num, "Error allocating memory for FFT data.\n");
		stop_reason = STOP_OUT_OF_MEM;
		goto exit;
	}
	cs.alt_x = gwalloc (&cs.gwdata);
	cs.u0 = gwalloc (&cs.gwdata);
	cs.d = gwalloc (&cs.gwdata);
	if (cs.alt_x == NULL || cs.u0 == NULL || cs.d == NULL) {
		OutputStr (thread_num, "Error allocating memory for Gerbicz error checking.\n");
		stop_reason = STOP_OUT_OF_MEM;
		goto exit;
	}

/* Format the string representation of the test number */

	strcpy (string_rep, gwmodulo_as_string (&cs.gwdata));

/* Init the title */

	sprintf (buf, "Cert %s", string_rep);
	title (thread_num, buf);

/* Loop reading from save files (and backup save files).  Limit number of backup */
/* files we try to read in case there is an error deleting bad save files. */

	readSaveFileStateInit (&read_save_file_state, thread_num, filename);
	for ( ; ; ) {

/* If there are no more save files, start off with the 1st cert squaring */

		if (! saveFileExists (&read_save_file_state)) {
			/* No save files existed, start from scratch. */
			cs.counter = 0;
			cs.error_count = 0;
			first_iter_msg = FALSE;
			break;
		}

/* Read a cert save file.  If successful, then if we've rolled back far enough for a restart break out of loop. */

		if (readCertSaveFile (&cs, read_save_file_state.current_filename, w)) {
			if (restart_counter < 0 || cs.counter <= (unsigned long) restart_counter) {
				first_iter_msg = TRUE;
				break;
			} else {
				// Don't treat save files we are ignoring due to a restart as bad.
				read_save_file_state.a_non_bad_save_file_existed = 0;
				_unlink (read_save_file_state.current_filename);
			}
		}

/* On read error, output message and loop to try the next backup save file. */

		else
			saveFileBad (&read_save_file_state);
	}

/* If this is a restart from an error, use the incremented error_count in restart_error_count */
/* rather than the error_count from a save file. */

	if (restart_error_count) cs.error_count = restart_error_count;

/* Output a message saying we are starting/resuming the certification. */
/* Also output the FFT length. */

	gwfft_description (&cs.gwdata, fft_desc);
	strcpy (buf, (cs.counter == 0) ? "Starting " : "Resuming ");
	sprintf (buf+strlen(buf), "Gerbicz error-checking certification of %s using %s\n", string_rep, fft_desc);
	OutputStr (thread_num, buf);

/* Clear all timers */

	clear_timers (timers, sizeof (timers) / sizeof (timers[0]));

/* Init vars for Test/Status and CommunicateWithServer */

	strcpy (w->stage, "CERT");
	inverse_explen = 1.0 / (double) w->cert_squarings;
	w->pct_complete = (double) cs.counter * inverse_explen;
	calc_interval_adjustments (&cs.gwdata, &output_adjustment, &title_adjustment);

/* If we are near the maximum exponent this fft length can test, then we */
/* will error check all iterations */

	near_fft_limit = exponent_near_fft_limit (&cs.gwdata);

/* Figure out the maximum round-off error we will allow.  By default this is 27/64 when near the FFT limit and 26/64 otherwise. */
/* We've found that this default catches errors without raising too many spurious error messages.  We let the user override */
/* this default for user "Never Odd Or Even" who tests exponents well beyond an FFT's limit.  He does his error checking by */
/* running the first-test and double-check simultaneously. */

	allowable_maxerr = IniGetFloat (INI_FILE, "MaxRoundoffError", (float) (near_fft_limit ? 0.421875 : 0.40625));

/* Get setting for verbosity of hardware error messages.  Force output of "confidence is excellent" when error checking. */

	error_count_messages = IniGetInt (INI_FILE, "ErrorCountMessages", 3);
	error_count_messages |= 0x8000;

/* Set the proper starting value and state if no save file was present */

	if (cs.counter == 0) {
		int	i, residue_size;
		uint32_t *array;		/* Array to contain the binary value */
		uint32_t arraylen;		/* Size of the array */
		char	md5[33], residue_md5[33];

		// Allocate an array for binary value
		residue_size = divide_rounding_up ((int) ceil(cs.gwdata.bit_length), 8);
		arraylen = divide_rounding_up (residue_size, 4);
		array = (uint32_t *) malloc (arraylen * sizeof(uint32_t));
		if (array == NULL) {
			OutputBoth (thread_num, "Error allocating memory for CERT starting value.\n");
			goto retry_work;
		}
		array[arraylen-1] = 0;		// Zero-pad the top few bytes

/* Contact PrimeNet for the starting residue */

		if (ProofGetData (w->assignment_uid, array, residue_size, md5) != PRIMENET_NO_ERROR || !isHex (md5)) {
			OutputBoth (thread_num, "Error getting CERT starting value.\n");
			free (array);
			goto retry_work;
		}
		md5_hexdigest_buffer (residue_md5, array, residue_size);
		strupper (md5);
		strupper (residue_md5);
		if (strcmp (md5, residue_md5) != 0) {
			OutputBoth (thread_num, "MD5 of downloaded starting value does not match.\n");
			free (array);
			goto retry_work;
		}

// BUG - find a way to commonize this code?? -- used in verifier, server_proof, proof_hash, read/write residue in prime95, etc

		// Convert from an array of bytes (LSB to MSB) to an array of uint32_t
		for (i = 0; i < (int) arraylen; i++) {
			uint32_t val;
			val = ((unsigned char *)&array[i])[3];
			val = (val << 8) + ((unsigned char *)&array[i])[2];
			val = (val << 8) + ((unsigned char *)&array[i])[1];
			val = (val << 8) + ((unsigned char *)&array[i])[0];
			array[i] = val;
		}

		// Convert binary to gwnum
		binarytogw (&cs.gwdata, array, arraylen, cs.x);
		free (array);

		// Generate a random initial shift count.  Initial shift count can't be larger than n.
		// Perform the initial rotate left.
		srand ((unsigned) time (NULL));
		cs.units_bit = ((rand () << 16) + rand ()) % w->n;
		gwrotate_left (&cs.gwdata, cs.x, cs.units_bit);

		// Gerbicz error checking expects alt_x to be a copy of x at the start of a block
		gwcopy (&cs.gwdata, cs.x, cs.alt_x);
		cs.alt_units_bit = cs.units_bit;
		cs.state = PRP_STATE_GERB_START_BLOCK;

		// Write a save file so that rolling back to the start does not require downloading the starting value again
		if (writeCertSaveFile (&cs, &write_save_file_state, w)) setWriteSaveFileSpecial (&write_save_file_state);
	}

/* Do the certification.  The certification is a long string of squarings, which is perfect for Gerbicz error checking. */
/* See the description of Gerbicz error checking in the prp routine.  In short, the x values we square are u[i] where */
/*	u[i]=u[0]^(2^i) mod mp */
/* The checksum function updated every L-th squaring is */
/*	d[t]=prod(i=0,t,u[L*i]) mod mp */
/* which can be calculated two different ways and compared to detect an error: */
/*	d[t]=d[t-1]*u[L*t] mod mp		(checksum #1) */
/*	d[t]=u[0]*d[t-1]^(2^L) mod mp		(checksum #2) */
/* The last few squarings that are too few to make a Gerbicz block are double-checked. */

	iters = 0;
	while (cs.counter < (unsigned long) w->cert_squarings) {
		gwnum	x;			/* Pointer to number to square */
		unsigned long *units_bit;	/* Pointer to units_bit to update */
		int	saving, saving_highly_reliable, actual_frequency;

/* If this is the first iteration of a Gerbicz error-checking block, then */
/* determine "L" -- the number of squarings between each Gerbicz multiplication */
/* We end this Gerbicz block after L^2 iterations.  */
/* If there aren't many iterations left, revert to simple double-checking. */

		if (cs.state == PRP_STATE_GERB_START_BLOCK) {
			int	iters_left = w->cert_squarings - cs.counter;
			if (iters_left < 49) {
				cs.state = PRP_STATE_DCHK_PASS1;
				cs.start_counter = cs.counter;
				cs.end_counter = cs.counter + iters_left;
				// Use a different shift count for alt_x so that the two double-check computations use different FFT data
				gwadd3o (&cs.gwdata, cs.alt_x, cs.alt_x, cs.alt_x, GWADD_FORCE_NORMALIZE); // Don't corrupt min/max roundoff error output
				cs.alt_units_bit = cs.alt_units_bit + 1;
				if (cs.alt_units_bit >= w->n) cs.alt_units_bit -= w->n;
			} else {
				int	gerbicz_block_size;
				double	adjustment;
				cs.state = PRP_STATE_GERB_MID_BLOCK;
				adjustment = IniGetFloat (INI_FILE, "PRPGerbiczCompareIntervalAdj", 1.0);
				if (adjustment < 0.001 || adjustment > 1.0) adjustment = 0.5;
				gerbicz_block_size = (int) (adjustment * IniGetInt (INI_FILE, "CertGerbiczCompareInterval", 100000));
				if (gerbicz_block_size < 25) gerbicz_block_size = 25;
				if (gerbicz_block_size > iters_left) gerbicz_block_size = iters_left;
				cs.L = (unsigned long) sqrt ((double) gerbicz_block_size);
				cs.start_counter = cs.counter;
				cs.next_mul_counter = cs.counter + cs.L;
				cs.end_counter = cs.counter + cs.L * cs.L;
				gwswap (cs.alt_x, cs.u0);		// Set u0 to a copy of x
				gwcopy (&cs.gwdata, cs.x, cs.d);	// Set d[0] to a copy of x
				if (IniGetInt (INI_FILE, "GerbiczVerbosity", 1) > 1) {
					sprintf (buf, "Start Gerbicz block of size %ld at iteration %ld.\n", cs.L * cs.L, cs.start_counter+1);
					OutputBoth (thread_num, buf);
				}
			}
		}

/* Save if we are stopping, right after we pass an errored iteration, several iterations before retesting */
/* an errored iteration so that we don't have to backtrack very far to do a gwsquare_carefully iteration */
/* (we don't do the iteration immediately before because a save operation may change the FFT data and make */
/* the error non-reproducible), and finally save if the save file timer has gone off. */

		stop_reason = stopCheck (thread_num);
		saving = stop_reason || cs.counter == last_counter-8 || cs.counter == last_counter || testSaveFilesFlag (thread_num);
		saving_highly_reliable = FALSE;

/* Round off error check before writing a save file, near an FFT size's limit, */
/* or check every iteration option is set, and every 128th iteration. */

		echk = ERRCHK || saving || near_fft_limit || (cs.counter & 127) == 0;
		gw_clear_maxerr (&cs.gwdata);

/* Do one certification squaring or Gerbicz multiply */

		timers[1] = 0.0;
		start_timer (timers, 1);

/* If we are doing one of the Gerbicz multiplies (not a squaring), then handle that here */

		if (cs.state == PRP_STATE_GERB_MID_BLOCK_MULT) {
			gwstartnextfft (&cs.gwdata, 0);		/* Do not start next forward FFT */
			gwsetnormroutine (&cs.gwdata, 0, 1, 0);	/* Always roundoff error check multiplies */
			gwsafemul (&cs.gwdata, cs.x, cs.d);	/* "Safe" multiply that does not change cs.x */
		} else if (cs.state == PRP_STATE_GERB_END_BLOCK_MULT) {
			gwstartnextfft (&cs.gwdata, 0);		/* Do not start next forward FFT */
			gwsetnormroutine (&cs.gwdata, 0, 1, 0);	/* Always roundoff error check multiplies */
			gwmul (&cs.gwdata, cs.u0, cs.alt_x);	/* Multiply to calc checksum #2.  u0 value can be destroyed. */
		} else if (cs.state == PRP_STATE_GERB_FINAL_MULT) {
			gwcopy (&cs.gwdata, cs.x, cs.u0);	// Copy x (before using it) for next Gerbicz block
			gwstartnextfft (&cs.gwdata, 0);		/* Do not start next forward FFT */
			gwsetnormroutine (&cs.gwdata, 0, 1, 0);	/* Always roundoff error check multiplies */
			gwsafemul (&cs.gwdata, cs.u0, cs.d);	/* "Safe" multiply to compute final d[t] value (checksum #1) */
		}

/* Otherwise, do a squaring iteration.  Use state to decide which number we are squaring. */

		else {
			if (cs.state == PRP_STATE_DCHK_PASS1 || cs.state == PRP_STATE_GERB_MID_BLOCK) {
				x = cs.x;
				units_bit = &cs.units_bit;
			} else {			// (cs.state == PRP_STATE_DCHK_PASS2 || cs.state == PRP_STATE_GERB_END_BLOCK) {
				x = cs.alt_x;
				units_bit = &cs.alt_units_bit;
			}

/* Decide if we can start the next forward FFT.  This is faster, but leaves the result in an "unsavable-to-disk" state. */

			gwstartnextfft (&cs.gwdata,
					!saving && !maxerr_recovery_mode && cs.counter != cs.next_mul_counter-1 &&
					cs.counter != cs.end_counter-1 && cs.counter != w->cert_squarings-1);
			gwsetnormroutine (&cs.gwdata, 0, echk, 0);

/* Square carefully during an error recovery.  This will protect us from roundoff errors up to (1.0 - 0.40625). */

			if (maxerr_recovery_mode && cs.counter == last_counter) {
				gwsquare_carefully (&cs.gwdata, x);
				maxerr_recovery_mode = 0;
				last_counter = 0xFFFFFFFF;
				echk = 0;
			} else
				gwsquare (&cs.gwdata, x);

			*units_bit <<= 1;
			if (*units_bit >= w->n) *units_bit -= w->n;
		}

// introduce an error every random # iterations when debugging highly reliable error checking
//#define INTRODUCE_ERRORS
#ifdef INTRODUCE_ERRORS
		if ((rand () & 0x7FFF) == 134)  // one out of 32768
			*cs.x += 5.0;
#endif

/* End iteration timing and increase count of iterations completed */

		end_timer (timers, 1);
		timers[0] += timers[1];
		iters++;

/* Update min/max round-off error */

		if (echk) {
			if (gw_get_maxerr (&cs.gwdata) < reallyminerr) reallyminerr = gw_get_maxerr (&cs.gwdata);
			if (gw_get_maxerr (&cs.gwdata) > reallymaxerr) reallymaxerr = gw_get_maxerr (&cs.gwdata);
		}

/* Check for excessive roundoff error.  Gerbicz error-checking will catch any problems caused by modest roundoff */
/* errors later.  However, if the round off error is really large, then results are certainly corrupt and we roll */
/* back immmediately.  If the large roundoff error is repeatable, then repeat the iteration using a safer, slower */
/* method.  This can happen when operating near the limit of an FFT. */

		if (echk && gw_get_maxerr (&cs.gwdata) > allowable_maxerr) {
			if (cs.counter == last_counter && gw_get_maxerr (&cs.gwdata) == last_maxerr) {
				OutputBoth (thread_num, ERROK);
				inc_error_count (3, &cs.error_count);
				gw_clear_error (&cs.gwdata);
				OutputBoth (thread_num, ERRMSG5);
				maxerr_recovery_mode = 1;
				restart_counter = cs.counter;		/* rollback to this iteration or earlier */
				sleep5 = FALSE;
				goto restart;
			} else {
				char	msg[100];
				sprintf (msg, ERRMSG1C, gw_get_maxerr (&cs.gwdata), allowable_maxerr);
				sprintf (buf, ERRMSG0, cs.counter+1, (long) w->cert_squarings, msg);
				OutputBoth (thread_num, buf);
				inc_error_count (1, &cs.error_count);
				if (gw_get_maxerr (&cs.gwdata) > IniGetFloat (INI_FILE, "RoundoffRollbackError", (float) 0.475)) {
					last_counter = cs.counter;
					last_maxerr = gw_get_maxerr (&cs.gwdata);
					restart_counter = cs.counter;		/* rollback to this iteration or earlier */
					sleep5 = FALSE;
					goto restart;
				}
			}
		}

/* Update counter, percentage complete */

		cs.counter++;
		w->pct_complete = (double) cs.counter * inverse_explen;

/* Output the title every so often */

		actual_frequency = (int) (ITER_OUTPUT * title_adjustment);
		if (actual_frequency < 1) actual_frequency = 1;
		if (cs.counter % actual_frequency == 0 || first_iter_msg) {
			sprintf (buf, "%.*f%% of cert %s", (int) PRECISION, trunc_percent (w->pct_complete), string_rep);
			title (thread_num, buf);
		}

/* Print a message every so often */

		actual_frequency = (int) (ITER_OUTPUT * output_adjustment);
		if (actual_frequency < 1) actual_frequency = 1;
		if ((cs.counter % actual_frequency == 0 && cs.state != PRP_STATE_GERB_MID_BLOCK_MULT &&
		     cs.state != PRP_STATE_GERB_END_BLOCK && cs.state != PRP_STATE_GERB_END_BLOCK_MULT &&
		     cs.state != PRP_STATE_GERB_FINAL_MULT) || first_iter_msg) {
			sprintf (buf, "Iteration: %ld / %ld [%.*f%%]",
				 cs.counter, (long) w->cert_squarings, (int) PRECISION, trunc_percent (w->pct_complete));
			/* Append a short form total errors message */
			if ((error_count_messages & 0xFF) == 1)
				make_error_count_message (cs.error_count, error_count_messages, buf + strlen (buf),
							  (int) (sizeof (buf) - strlen (buf)));
			/* Truncate first message */
			if (first_iter_msg) {
				strcat (buf, ".\n");
				clear_timer (timers, 0);
				first_iter_msg = FALSE;
			}
			/* In v28.5 and later, format a consise message including the ETA */
			else if (!CLASSIC_OUTPUT) {
				double speed;
				/* Append roundoff error */
				if ((OUTPUT_ROUNDOFF || ERRCHK) && reallymaxerr >= 0.001) {
					sprintf (buf+strlen(buf), ", roundoff: %5.3f", reallymaxerr);
					if (!CUMULATIVE_ROUNDOFF) reallyminerr = 1.0, reallymaxerr = 0.0;
				}
				/* Append ms/iter */
				speed = timer_value (timers, 0) / (double) iters;
				sprintf (buf+strlen(buf), ", ms/iter: %6.3f", speed * 1000.0);
				clear_timer (timers, 0);
				iters = 0;
				/* Append ETA */
				formatETA ((w->cert_squarings - cs.counter) * speed, buf+strlen(buf));
				strcat (buf, "\n");
			}
			/* Format the classic (pre-v28.5) message */
			else {
				/* Append optional roundoff message */
				if (ERRCHK && cs.counter > 30) {
					sprintf (buf+strlen(buf), ".  Round off: %10.10f to %10.10f", reallyminerr, reallymaxerr);
					if (!CUMULATIVE_ROUNDOFF) reallyminerr = 1.0, reallymaxerr = 0.0;
				}
				if (CUMULATIVE_TIMING) {
					strcat (buf, ".  Total time: ");
					print_timer (timers, 0, buf, TIMER_NL);
				} else {
					strcat (buf, ".  Per iteration time: ");
					divide_timer (timers, 0, iters);
					print_timer (timers, 0, buf, TIMER_NL | TIMER_CLR);
					iters = 0;
				}
			}
			OutputStr (thread_num, buf);

/* Output a verbose message showing the error counts.  This way a user is likely to */
/* notice a problem without reading the results.txt file. */

			if ((error_count_messages & 0xFF) >= 2 &&
			    make_error_count_message (cs.error_count, error_count_messages, buf, sizeof (buf)))
				OutputStr (thread_num, buf);
		}

/* If double-checking, at end of pass 1 rollback counter and start computing alt_x. */
/* If double-checking, at end of pass 2 compare values and move onto next block. */

		if (cs.state == PRP_STATE_DCHK_PASS1) {
			if (cs.counter < cs.end_counter);		// Do next iteration
			else if (cs.counter == cs.end_counter) {	// Switch to alt_x computations
				cs.state = PRP_STATE_DCHK_PASS2;
				cs.counter = cs.start_counter;
			} else {					// Can't happen
				OutputBoth (thread_num, ERRMSG9);
				inc_error_count (6, &cs.error_count);
				restart_counter = -1;			/* rollback to any save file */
				sleep5 = FALSE;
				goto restart;
			}
		}
		if (cs.state == PRP_STATE_DCHK_PASS2) {
			if (cs.counter < cs.end_counter);		// Do next iteration
			else if (cs.counter == cs.end_counter) {	// Compare x and alt_x
				if (!areTwoPRPValsEqual (&cs.gwdata, w->n, cs.x, cs.units_bit, cs.alt_x, cs.alt_units_bit)) {
					sprintf (buf, ERRMSG6, cs.start_counter);
					OutputBoth (thread_num, buf);
					inc_error_count (7, &cs.error_count);
					restart_counter = cs.start_counter;		/* rollback to this iteration */
					sleep5 = FALSE;
					goto restart;
				}
				cs.state = PRP_STATE_GERB_START_BLOCK;
			} else {					// Can't happen
				OutputBoth (thread_num, ERRMSG9);
				inc_error_count (6, &cs.error_count);
				restart_counter = -1;			/* rollback to any save file */
				sleep5 = FALSE;
				goto restart;
			}
		}

/* Handle all the possible Gerbicz states.  See if this is an L-th iteration that needs to do a checksum multiply. */
/* Also check if this is the L^2-th iteration where we switch to an alternate method to compute checksum #2. */
/* NOTE: We do not need to worry about shift_counts in the checksum as both checksums will end up with the same shift count. */

		// Just did a normal cert squaring
		if (cs.state == PRP_STATE_GERB_MID_BLOCK) {
			if (cs.counter < cs.next_mul_counter);		// Do next iteration
			else if (cs.counter == cs.end_counter) {	// Delay last checksum #1 multiply, start checksum #2 calculation
				if (IniGetInt (INI_FILE, "GerbiczVerbosity", 1) > 1) OutputStr (thread_num, "Start Gerbicz error check.\n");
				// At end of Gerbicz block, switch to "L" squarings of alt_x to create Gerbicz checksum #2 value
				gwcopy (&cs.gwdata, cs.d, cs.alt_x);	// Copy d[t-1] to alt_x
				cs.state = PRP_STATE_GERB_END_BLOCK;	// Squaring alt_x state
				cs.counter -= cs.L;			// L squarings
			} else if (cs.counter == cs.next_mul_counter) {	// Do a checksum #1 multiply next
				// Back counter up by one and do one multiply in the computation of Gerbicz checksum #1 value
				cs.state = PRP_STATE_GERB_MID_BLOCK_MULT;
				cs.counter -= 1;
			} else {					// Can't happen
				OutputBoth (thread_num, ERRMSG9);
				inc_error_count (6, &cs.error_count);
				restart_counter = -1;			/* rollback to any save file */
				sleep5 = FALSE;
				goto restart;
			}
		}

		// Just did a a checksum #1 multiply at the end of a block of L normal squarings
		else if (cs.state == PRP_STATE_GERB_MID_BLOCK_MULT) {
			if (cs.counter < cs.end_counter) {		// In middle of Gerbicz block, do another "L" squarings
				cs.state = PRP_STATE_GERB_MID_BLOCK;
				cs.next_mul_counter += cs.L;
			} else {					// Can't happen
				OutputBoth (thread_num, ERRMSG9);
				inc_error_count (6, &cs.error_count);
				restart_counter = -1;			/* rollback to any save file */
				sleep5 = FALSE;
				goto restart;
			}
		}

		// Just did a checksum #2 squaring
		else if (cs.state == PRP_STATE_GERB_END_BLOCK) {
			if (cs.counter < cs.end_counter);		// Do next iteration in computing checksum #2
			else if (cs.counter == cs.end_counter) {	// Next do final multiply in computing checksum #2
				cs.state = PRP_STATE_GERB_END_BLOCK_MULT;
				cs.counter -= 1;
			} else {					// Can't happen
				OutputBoth (thread_num, ERRMSG9);
				inc_error_count (6, &cs.error_count);
				restart_counter = -1;			/* rollback to any save file */
				sleep5 = FALSE;
				goto restart;
			}
		}

		// Just did the checksum #2 multiply at the end of L checksum #2 squarings
		else if (cs.state == PRP_STATE_GERB_END_BLOCK_MULT) {
			if (cs.counter == cs.end_counter) {		// Next do final multiply in computing checksum #1
				cs.state = PRP_STATE_GERB_FINAL_MULT;
				cs.counter -= 1;
			} else {					// Can't happen
				OutputBoth (thread_num, ERRMSG9);
				inc_error_count (6, &cs.error_count);
				restart_counter = -1;			/* rollback to any save file */
				sleep5 = FALSE;
				goto restart;
			}
		}

		// Just did the final checksum #1 multiply which we delayed until checksum #2 completed
		else if (cs.state == PRP_STATE_GERB_FINAL_MULT) {
			double	gerbicz_block_size_adjustment;
			// We adjust the compare interval size downward when errors occur and upwards when they dont.
			// That way buggy machines will lose fewer iterations when rolling back.
			gerbicz_block_size_adjustment = IniGetFloat (INI_FILE, "PRPGerbiczCompareIntervalAdj", 1.0);
			if (gerbicz_block_size_adjustment < 0.001 || gerbicz_block_size_adjustment > 1.0) gerbicz_block_size_adjustment = 0.5;
			// Compare alt_x, d (the two Gerbicz checksum values that must match)
			if (!areTwoPRPValsEqual (&cs.gwdata, w->n, cs.alt_x, 0, cs.d, 0)) {
				sprintf (buf, ERRMSG7, cs.start_counter);
				OutputBoth (thread_num, buf);
				gerbicz_block_size_adjustment *= 0.25;		/* This will halve next L */
				if (gerbicz_block_size_adjustment < 0.001) gerbicz_block_size_adjustment = 0.001;
				IniWriteFloat (INI_FILE, "PRPGerbiczCompareIntervalAdj", (float) gerbicz_block_size_adjustment);
				inc_error_count (7, &cs.error_count);
				restart_counter = cs.start_counter;		/* rollback to this iteration */
				sleep5 = FALSE;
				goto restart;
			}
			if (IniGetInt (INI_FILE, "GerbiczVerbosity", 1)) {
				sprintf (buf, "Gerbicz error check passed at iteration %ld.\n", cs.counter);
				OutputStr (thread_num, buf);
			}
			gerbicz_block_size_adjustment *= 1.0473;		/* 30 good blocks to double L */
			if (gerbicz_block_size_adjustment > 1.0) gerbicz_block_size_adjustment = 1.0;
			IniWriteFloat (INI_FILE, "PRPGerbiczCompareIntervalAdj", (float) gerbicz_block_size_adjustment);
			/* Start next Gerbicz block.  Both x and alt_x must be identical at start of next block. */
			cs.state = PRP_STATE_GERB_START_BLOCK;
			gwswap (cs.alt_x, cs.u0);
			cs.alt_units_bit = cs.units_bit;
			/* We've reached a verified iteration, create a save file and mark it highly reliable. */
			/* But if there are less than 1000 iterations left on a reliable machine, don't bother creating the save file. */
			if (w->cert_squarings - cs.counter >= 1000 || gerbicz_block_size_adjustment < 1.0) {
				saving = TRUE;
				saving_highly_reliable = TRUE;
			}
		}

/* Write a save file every DISK_WRITE_TIME minutes */

		if (saving) {
			if (writeCertSaveFile (&cs, &write_save_file_state, w)) {
				// Mark save files that contain verified computations.  This will keep the save file
				// for a longer period of time (i.e. will not be replaced by a save file that does
				// not also contain verified computations).
				if (saving_highly_reliable) setWriteSaveFileSpecial (&write_save_file_state);
			}
		}

/* If an escape key was hit, write out the results and return */

		if (stop_reason) {
			sprintf (buf, "Stopping certification of %s at iteration %ld [%.*f%%]\n",
				 string_rep, cs.counter, (int) PRECISION, trunc_percent (w->pct_complete));
			OutputStr (thread_num, buf);
			goto exit;
		}
	}

/* Make sure cert state is valid.  We cannot be in the middle of a double-check or in the middle of a Gerbicz block */

	if (cs.state != PRP_STATE_GERB_START_BLOCK) {
		OutputBoth (thread_num, ERRMSG9);
		inc_error_count (6, &cs.error_count);
		restart_counter = -1;			/* rollback to any save file */
		sleep5 = FALSE;
		goto restart;
	}

/* Undo the shift count, generate the 256-bit SHA-3 hash */

	if (!gwrotate_right (&cs.gwdata, cs.x, cs.units_bit) || !sha3_gwnum (&cs.gwdata, cs.x, &hash)) {
		OutputBoth (thread_num, ERRMSG8);
		inc_error_count (2, &cs.error_count);
		restart_counter = -1;			/* rollback to any save file */
		sleep5 = TRUE;
		goto restart;
	}

/* Free up some memory */

	gwfree (&cs.gwdata, cs.x);
	gwfree (&cs.gwdata, cs.alt_x);
	gwfree (&cs.gwdata, cs.u0);
	gwfree (&cs.gwdata, cs.d);

/* Print results */

	sprintf (buf, "%s certification hash value %s.", string_rep, hash_to_string (hash));
	sprintf (buf+strlen(buf), " Wh%d: %08lX,", PORT, SEC1 (w->n));
	if (cs.units_bit) sprintf (buf+strlen(buf), "%ld,", cs.units_bit);
	sprintf (buf+strlen(buf), "%08lX\n", cs.error_count);
	OutputStr (thread_num, buf);
	formatMsgForResultsFile (buf, w);

/* Update the output file */

	if (IniGetInt (INI_FILE, "OutputCerts", 0))
		writeResults (buf);

/* Format a JSON version of the result.  An example follows: */
/* {"exponent":25000000, "worktype":"Cert", "sha3-hash":"0123456789ABCDEF0123456789ABCDEF0123456789ABCDEF0123456789ABCDEF", */
/* "fft-length":4096000, "squarings":100000, "error-code":"00010000", */
/* "security-code":"C6B0B26C", "program":{"name":"prime95", "version":"30.1", "build":"1"}, "timestamp":"2019-01-15 23:28:16", */
/* "errors":{"gerbicz":0}, "user":"gw_2", "cpu":"basement", "aid":"FF00AA00FF00AA00FF00AA00FF00AA00"} */

	sprintf (JSONbuf, "{");
	sprintf (JSONbuf+strlen(JSONbuf), "\"worktype\":\"Cert\"");
	JSONaddExponent (JSONbuf, w);
	sprintf (JSONbuf+strlen(JSONbuf), ", \"sha3-hash\":\"%s\"", hash_to_string (hash));
	sprintf (JSONbuf+strlen(JSONbuf), ", \"fft-length\":%lu", cs.gwdata.FFTLEN);
	sprintf (JSONbuf+strlen(JSONbuf), ", \"squarings\":%d", w->cert_squarings);
	if (cs.units_bit) sprintf (JSONbuf+strlen(JSONbuf), ", \"shift-count\":%ld", cs.units_bit);
	sprintf (JSONbuf+strlen(JSONbuf), ", \"error-code\":\"%08lX\"", cs.error_count);
	sprintf (JSONbuf+strlen(JSONbuf), ", \"security-code\":\"%08lX\"", SEC1(w->n));
	JSONaddProgramTimestamp (JSONbuf);
	sprintf (JSONbuf+strlen(JSONbuf), ", \"errors\":{\"gerbicz\":%lu}", (cs.error_count >> 20) & 0xF);
	JSONaddExponentKnownFactors (JSONbuf, w);
	JSONaddUserComputerAID (JSONbuf, w);
	strcat (JSONbuf, "}");
	if (IniGetInt (INI_FILE, "OutputJSONCerts", 1)) writeResultsJSON (JSONbuf);

/* Output results to the server */

	{
		struct primenetAssignmentResult pkt;
		memset (&pkt, 0, sizeof (pkt));
		strcpy (pkt.computer_guid, COMPUTER_GUID);
		strcpy (pkt.assignment_uid, w->assignment_uid);
		strcpy (pkt.message, buf);
		pkt.result_type = PRIMENET_AR_CERT;
		pkt.k = w->k;
		pkt.b = w->b;
		pkt.n = w->n;
		pkt.c = w->c;
		strcpy (pkt.cert_hash, hash_to_string (hash));
		sprintf (pkt.error_count, "%08lX", cs.error_count);
		pkt.shift_count = cs.units_bit;
		pkt.fftlen = gwfftlen (&cs.gwdata);
		pkt.done = TRUE;
		strcpy (pkt.JSONmessage, JSONbuf);
		spoolMessage (PRIMENET_ASSIGNMENT_RESULT, &pkt);
	}

/* Delete the save files and report work is done */

	unlinkSaveFiles (&write_save_file_state);
	goto work_complete;

/* Retry up to 8 times.  Output retry work later message and return. */

retry_work:
	echk = IniSectionGetInt (INI_FILE, SEC_Internals, KEY_CertErrorCount, 0) + 1;
	IniSectionWriteInt (INI_FILE, SEC_Internals, KEY_CertErrorCount, echk);
	if (echk >= 8) goto abandon_work;
	OutputStr (thread_num, "Will retry certification later.\n");
	stop_reason = STOP_RETRY_LATER;
	goto exit;

/* Abandon work and exit */

abandon_work:
	sprintf (buf, "Abandoning certification of M%lu.\n", w->n);
	OutputBoth (thread_num, buf);
	unreserve (w->n);
	goto work_complete;

/* Return work unit complete "error" code */

work_complete:
	IniSectionWriteString (INI_FILE, SEC_Internals, KEY_CertErrorCount, NULL);
	stop_reason = STOP_WORK_UNIT_COMPLETE;
	goto exit;

/* Cleanup and exit */

exit:	gwdone (&cs.gwdata);
	return (stop_reason);

/* An error occurred, output a message saying we are restarting, sleep, */
/* then try restarting at last save point. */

restart:if (sleep5) OutputBoth (thread_num, ERRMSG2);
	OutputBoth (thread_num, ERRMSG3);

/* Save the incremented error count to be used in the restart rather than the error count read from a save file */

	restart_error_count = cs.error_count;

/* Sleep five minutes before restarting */

	if (sleep5) {
		stop_reason = SleepFive (thread_num);
		if (stop_reason) {
			gwdone (&cs.gwdata);
			return (stop_reason);
		}
	}

/* Return so that last continuation file is read in */

	gwdone (&cs.gwdata);
	goto begin;
}
