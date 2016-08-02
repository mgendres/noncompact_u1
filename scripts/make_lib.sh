#!/bin/bash

LIB=noncompact_u1
SRCDIR=$PWD
BUILDDIR=$PWD
SRCDIR2=$SRCDIR/src
BUILDDIR2=$BUILDDIR/objs
FFTWDIR=$HOME/Documents/libraries/fftw-3.3.3

CC=gcc
CXX=g++
AR=ar

LDFLAGS=""
CXXFLAGS="-O2"
ARFLAGS="ruv"

INCLUDE_FLAGS="-I$SRCDIR/include -I$FFTWDIR/include"

if [ -d $BUILDDIR2 ]; then
  rm -r $BUILDDIR2
fi

mkdir $BUILDDIR2
cd $BUILDDIR2

$CXX $CXXFLAGS $LDFLAGS $INCLUDE_FLAGS -c $SRCDIR2/*.C
$AR $ARFLAGS $BUILDDIR/$LIB.a $BUILDDIR2/*.o
ranlib  $BUILDDIR/$LIB.a
