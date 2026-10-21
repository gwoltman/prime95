#!/bin/sh
# Build the official Linux release of mprime, on the build machine itself.
#
# Produces, in the output directory (default: prime95-release-<version> next
# to the source tree):
#   p95v<version>.linux64.tar.gz   on 64-bit Linux
#   p95v<version>.linux32.tar.gz   on 32-bit Linux
#   SHA256SUMS                     (added to, if it already exists)
# where <version> comes from commonc.h, for example 3107b02.  Each tarball
# holds mprime, the GMP shared library it was linked against and the text
# files.
#
# It builds the source tree it is in, which should be fresh: a new clone, or
# what release/vm.ps1 sends.  It adds the private security files, builds
# gwnum, factor and mprime, checks that mprime is a trusted build of the right
# version, and packages it.  The security files are deleted from the tree
# when it exits.
#
# Needs everything BUILDING.md lists for Linux; uasm and objconv must be on
# PATH.
#
#   release/unix.sh [-s security-dir] [-o output-dir] [-b boost-dir] [-g gmp-lib-dir]
#
#   -s  where security.h, security.c and secure5.c are
#       (default: prime95-security next to the source tree)
#   -b  the Boost headers (default: ~/boost_1_76_0)
#   -g  the directory holding libgmp.so (default: /usr/local/lib)

set -e

tree=$(cd "$(dirname "$0")/.." && pwd)
secdir=$tree/../prime95-security
outdir=
boost=$HOME/boost_1_76_0
gmplib=/usr/local/lib
while getopts s:o:b:g: opt; do
	case $opt in
	s) secdir=$OPTARG ;;
	o) outdir=$OPTARG ;;
	b) boost=$OPTARG ;;
	g) gmplib=$OPTARG ;;
	*) exit 2 ;;
	esac
done

step() { echo "== $*"; }
die() { echo "release/unix.sh: $*" >&2; exit 1; }

case $(uname -s)-$(uname -m) in
Linux-x86_64)	platform=linux64; appname=Linux64; asm=makeasm64; gwmake=make64; mprimemake=make64 ;;
Linux-i?86)	platform=linux32; appname=Linux; asm=makeasm32; gwmake=makefile; mprimemake=makefile ;;
*)		die "no release build is defined for $(uname -s) $(uname -m)" ;;
esac
command -v uasm > /dev/null || die "uasm is not on PATH"
command -v objconv > /dev/null || die "objconv is not on PATH"
[ -f "$boost/boost/version.hpp" ] || die "Boost headers not found in $boost (use -b)"
[ -f "$gmplib/libgmp.so" ] || die "libgmp.so not found in $gmplib (use -g)"

# Security files
secfiles="security.h security.c secure5.c"
for f in $secfiles; do
	[ -f "$secdir/$f" ] || die "$f not found in $secdir (use -s)"
done
grep -q SECURITY_MODULES_PRESENT "$secdir/security.h" || die "security.h in $secdir does not define SECURITY_MODULES_PRESENT"

# Version, from commonc.h
major=$(sed -n 's/^#define[ 	]*VERSION[ 	]*"\([0-9]*\)\.[0-9]*".*/\1/p' "$tree/commonc.h")
minor=$(sed -n 's/^#define[ 	]*VERSION[ 	]*"[0-9]*\.\([0-9]*\)".*/\1/p' "$tree/commonc.h")
build=$(sed -n 's/^#define[ 	]*BUILD_NUM[ 	]*"\([0-9]*\)".*/\1/p' "$tree/commonc.h")
[ -n "$major" ] && [ -n "$minor" ] && [ -n "$build" ] || die "VERSION or BUILD_NUM not found in commonc.h"
ver=$(printf '%s%02db%02d' "$major" "$minor" "$build")
[ -n "$outdir" ] || outdir=$tree/../prime95-release-$ver
echo "Version $major.$minor build $build ($ver), $platform"

[ ! -e "$tree/linux/mprime" ] || die "$tree is not a fresh tree: linux/mprime already exists"
cleanup() { for f in $secfiles; do rm -f "$tree/$f"; done; }
trap cleanup EXIT
trap 'exit 1' HUP INT TERM
for f in $secfiles; do cp "$secdir/$f" "$tree/"; done

# Build.  Compiler output goes to release/build-<platform>.log.
jobs=$(getconf _NPROCESSORS_ONLN 2>/dev/null || echo 2)
log=$tree/release/build-$platform.log
: > "$log"
run() {
	dir=$1; shift
	step "$dir: make $*"
	(cd "$tree/$dir" && make "$@") >> "$log" 2>&1 || { tail -30 "$log"; die "make $* in $dir failed (see $log)"; }
}
run gwnum -j"$jobs" -f $asm
run gwnum -j"$jobs" -f $gwmake
run prime95 -f $asm
run linux -j"$jobs" -f $mprimemake BOOST="$boost"

# Check that mprime is a trusted build of this version
mprime=$tree/linux/mprime
! grep -q 'Untrusted Prime95' "$mprime" || die "mprime is an untrusted build: the security files were not used"
app=$("$mprime" -v)
echo "$app"
case $app in
*": $appname,Prime95,v$major.$minor,build $build") ;;
*) die "mprime -v should report $appname,Prime95,v$major.$minor,build $build" ;;
esac

# Package
step "Packaging into $outdir"
mkdir -p "$outdir"
outdir=$(cd "$outdir" && pwd)
stage=$tree/release/stage-$platform
rm -rf "$stage"
mkdir "$stage"
cp "$mprime" "$stage/"
chmod 755 "$stage/mprime"
# The real library file and its symbolic links, as installed
real=$(cd "$gmplib" && ls -l libgmp.so | sed 's/.* -> //')
case $real in
libgmp.so.*) ;;
*) die "$gmplib/libgmp.so is not a link to libgmp.so.<version>" ;;
esac
cp "$gmplib/$real" "$stage/"
chmod 755 "$stage/$real"
ln -s "$real" "$stage/libgmp.so"
ln -s "$real" "$stage/libgmp.so.10"
# The text files.  prime95's have always shipped with CRLF line endings;
# linux/whatsnew.txt is the Linux list of changes.
for f in license.txt readme.txt stress.txt undoc.txt; do
	awk '{ sub(/\r$/, ""); printf "%s\r\n", $0 }' "$tree/prime95/$f" > "$stage/$f"
done
cp "$tree/linux/whatsnew.txt" "$stage/"
chmod 644 "$stage"/*.txt
tarball=p95v$ver.$platform.tar.gz
(cd "$stage" && tar -czf "$outdir/$tarball" $(ls | sort))
rm -rf "$stage"
(cd "$outdir" && {
	[ ! -f SHA256SUMS ] || grep -v "  $tarball\$" SHA256SUMS > SHA256SUMS.new || true
	sha256sum "$tarball" >> SHA256SUMS.new
	sort -k 2 SHA256SUMS.new > SHA256SUMS
	rm SHA256SUMS.new
})
echo "  $tarball: $(tar -tzf "$outdir/$tarball" | wc -l) files, $(du -k "$outdir/$tarball" | cut -f1) KB"
echo "Done: $ver for $platform"
