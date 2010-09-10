def test_bargraph():
    """Playing around with bargraphs and facets"""

    with open('test_data', 'w') as f:
        f.write('Sample\tValue\n')
        for idx, val in enumerate('glioblastoma'):
            f.write('%d\t%d\n' % (idx*2, idx))

test_bargraph()
