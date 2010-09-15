"""Testing bed_tools.py"""
import sys, nose.tools
sys.path.append('../')
import bed_tools

def test_load_bed():
    """See if I get the right tracks"""

    bed_file = '../data/nimblegen/2.1M_Human_Exome_Annotation/2.1M_Human_Exome.bed'
    chr2st2end, chr2posLs = bed_tools.load_bed(bed_file, 'Target Regions')
    chr_count = len(chr2st2end.keys())
    nose.tools.assert_equal(chr_count, 24,
                            'dont have 24 chrs. have ' + str(chr_count))
    for chr in chr2posLs:
        for idx in xrange(len(chr2posLs[chr])-1):
            nose.tools.assert_true(chr2posLs[chr][idx] < 
                                   chr2posLs[chr][idx+1],
                                   'chr2posLs is not sorted')
        nose.tools.assert_true(len(chr2posLs[chr]) > 400,
                               chr + ' length is ' + str(len(chr2posLs[chr])))

def test_find_location_in_bed():
    """See if I'm using searching correctly w/ binary search"""

    # this is in the file
    chr = 'chr5'
    pos = 180603265
    bed_file = '../data/nimblegen/2.1M_Human_Exome_Annotation/2.1M_Human_Exome.bed'
    chr2st2end, chr2posLs = bed_tools.load_bed(bed_file, 'NimbleGen Tiled Regions')
    location = bed_tools.find_location_in_bed(chr, pos, chr2posLs, chr2st2end)
    nose.tools.assert_equal(location, 180603263,
                            'did not find location, found '
                            + str(location))

    # this is not in the file
    # too large
    chr = 'chrY'
    pos = 26180101
    location = bed_tools.find_location_in_bed(chr, pos, chr2posLs, chr2st2end)
    nose.tools.assert_true(not location,
                           'found location in too big test, but should not have, found '
                           + str(location))

    # this is not in the file
    # too small
    chr = 'chr12'
    pos = 100
    location = bed_tools.find_location_in_bed(chr, pos, chr2posLs, chr2st2end)
    nose.tools.assert_true(not location,
                           'found location in too small test, but should not have, found '
                           + str(location))

    # this is not in the file
    # between capture regions
    chr = 'chr6'
    pos = 170731150
    location = bed_tools.find_location_in_bed(chr, pos, chr2posLs, chr2st2end)
    nose.tools.assert_true(not location,
                           'found location in between test, but should not have, found '
                           + str(location))

def test_find_ge():
    """Can I use the bisect module correctly?"""

    ls = [1,2,3,4,5]
    val, idx = bed_tools.find_ge(ls,2)
    nose.tools.assert_equal(2, val, 'not using bisect right. '
                            + str(val) + ' was found')
    val, idx = bed_tools.find_ge(ls,0)
    nose.tools.assert_equal(val, 1, 'did not find start of list')
    val, idx = bed_tools.find_ge(ls,6)
    nose.tools.assert_equal(val, 5, 'did not find end of list')
