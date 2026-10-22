# Prime95 and mprime

Prime95 (Windows) and mprime (Linux, FreeBSD and macOS) are the client
software of the [Great Internet Mersenne Prime Search](https://www.mersenne.org)
(GIMPS).  They search for Mersenne primes, primes of the form 2^p - 1, by
running probable-prime tests, Lucas-Lehmer tests, trial factoring, P-1 and
ECM, and get their work from and report results to the PrimeNet server.  They
are also widely used to stress test computer hardware.

Most people should download a ready-built release from
https://www.mersenne.org/download/ instead of building from this repository.

The arithmetic is done by **gwnum**, a library of highly optimized FFT
multiplication of large numbers, written largely in x86 assembly.  Other
programs can build on it; see [The gwnum package](BUILDING.md#the-gwnum-package).

## Building

[BUILDING.md](BUILDING.md) explains how to build on Windows, Linux (64-bit and
32-bit), FreeBSD and macOS, and how the official releases are made.

Builds from this repository identify themselves to PrimeNet as "Untrusted
Prime95": the code that marks a build as official is kept private.  They are
otherwise complete.

## Layout

| Directory | Contents |
|---|---|
| (top level) | Code shared by Prime95 and mprime: work scheduling, PrimeNet communication, ECM, P-1, proofs |
| `gwnum/` | The gwnum library |
| `prime95/` | The Windows GUI, the trial factoring assembly code, and the text files shipped with every release |
| `linux/` | mprime, with makefiles for Linux, FreeBSD and macOS |
| `winnt/` | The Windows service version (ntprime) |
| `macosx/` | The old macOS GUI (not currently built) |
| `release/` | The scripts that build the official releases |
| `third_party/` | Code by others: cJSON, SQLite, QD, MD5 and SHA-3 |

## License

See [LICENSE](LICENSE).  Use of this software is governed by the GIMPS terms
at https://www.mersenne.org/legal/, including the rules for prize money for
new Mersenne primes.

Bug reports and suggestions are welcome at the
[mersenneforum.org](https://www.mersenneforum.org) software forum.
