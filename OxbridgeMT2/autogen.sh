#!/bin/sh -e
test -n "$srcdir" || srcdir=`dirname "$0"`
test -n "$srcdir" || srcdir=.

autoreconf --force --install --verbose "$srcdir"
rm INSTALL
svn update INSTALL

test -n "$NOCONFIGURE" || "$srcdir/configure" "$@"
