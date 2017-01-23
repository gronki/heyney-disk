#!/usr/bin/env python

from sys import argv
from pyminiconf import pyminiconf_read
from pickle import dump as pidump
from tarfile import open as taropen
from gzip import open as gzopen
from re import search,sub
from re import compile as re_compile

def readcd(fc,fd):
    from numpy import dtype,loadtxt
    dt = dtype([(s[:16].strip().lower(), s[16:].strip(), 1) for s in fc])
    return loadtxt(fd, skiprows=2, dtype=dt)

def col2python(fn):
    rexp_tgz = re_compile(r'(\.tar\.gz|\.tgz)$')
    rexp_dat = re_compile(r'\.dat$')
    if search(rexp_tgz,fn) != None:
        with taropen(fn,"r:gz") as tar:
            fc = tar.extractfile(sub(rexp_tgz,'.col',fn))
            fd = tar.extractfile(sub(rexp_tgz,'.dat',fn))
            ft = tar.extractfile(sub(rexp_tgz,'.txt',fn))
            return readcd(fc,fd), pyminiconf_read(ft)
    elif search(rexp_dat,fn) != None:
        with    open(fn,'r') as fd, \
                open(sub(rexp_dat,'.col',fn),'r') as fc, \
                open(sub(rexp_dat,'.txt',fn),'r') as ft:
            return readcd(fc,fd), pyminiconf_read(ft)
    else:
        raise Exception("It should be .tar.gz, .tgz or .dat file, got: %s." % fn)

if __name__ == '__main__':
    rexp_allowed = re_compile(r'(\.tar\.gz|\.tgz|\.dat)$')
    for fn in argv[1:]:
        try:
            p = col2python(fn)
            if p != None:
                pidump(p, gzopen(sub(rexp_allowed,'.pickle.gz',fn),'wb'))
        except Exception as e:
            print "col2python error: %s" % e
