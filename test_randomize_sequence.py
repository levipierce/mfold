from randomize_sequence import randomize_sequence as rs
import unittest


class TestRandomizeSequenceFunctions(unittest.TestCase):

    def setUp(self):
        self.seq_fname = "HSPA1B-mfoldinterval.fa"

    def test_read_sequence(self):
        expected_length = 2099
        # make sure the randomized sequence does not lose any elements
        seq, maps = rs.read_sequence(self.seq_fname)
        self.assertEqual(expected_length, len(seq))
        # relative to absolute map
        self.assertEqual(expected_length, len(maps[0].keys()))
        # absolute to relative map
        self.assertEqual(expected_length, len(maps[1].keys()))

    def test_no_masking(self):
        seq, maps = rs.read_sequence(self.seq_fname)
        rsr = rs.shuffle_sequence(seq, maps, n=100, mask=None)
        similarity_count = 0
        for idx, s in enumerate(rsr.seq):
            if seq[idx] == s:
                similarity_count+=1
        print "similarity count : ", similarity_count
        self.assertLess(similarity_count,len(seq))

    def test_masking_snp(self):
        """
        HSE1:		31795312-31795325
        HSE2:		31795405-31795418
        TATA:		31795486-31795491
        TSS: 		31795512-31795512
        TRIM28: 	31795571-31795578
        PB1: 		31795743-31795749
        PB2: 		31796358-31796364
        """
        mask={"snp":[(31795512,31795512)]}

        seq, maps = rs.read_sequence(self.seq_fname)
        abs_rel_map = maps[1]
       
        snp = abs_rel_map[31795512]
        print seq[snp]
        rsr = rs.shuffle_sequence(seq, maps, n=100, mask=("snp", mask["snp"])
            )
        print rsr.seq[snp]
        self.assertEqual(seq[snp], rsr.seq[snp])

    def test_masking_region(self):
        """
        HSE1:		31795312-31795325
        HSE2:		31795405-31795418
        TATA:		31795486-31795491
        TSS: 		31795512-31795512
        TRIM28: 	31795571-31795578
        PB1: 		31795743-31795749
        PB2: 		31796358-31796364
        """
        mask={"HSE1":[(31795312,31795325)]}

        seq, maps = rs.read_sequence(self.seq_fname)
        abs_rel_map = maps[1]
       
        expected_seq = ''.join([seq[abs_rel_map[s]] for s in range(31795312,31795325+1)])
        rsr = rs.shuffle_sequence(seq, maps, n=100, mask=("HSE1", mask["HSE1"])
            )
        region_seq =  ''.join([rsr.seq[abs_rel_map[s]] for s in range(31795312,31795325+1)])
        self.assertEqual(expected_seq, region_seq)

        #similarity count should be higher too...
    def test_masking_first_bp(self):
        mask = {"HSE1" : [(31795202, 31795202)]}

        seq, maps = rs.read_sequence(self.seq_fname)
        abs_rel_map = maps[1]

        expected_seq = ''.join([seq[abs_rel_map[s]] for s in range(31795202, 31795202+1)])
        print expected_seq
        rsr = rs.shuffle_sequence(seq, maps, n=100, mask=("HSE1", mask["HSE1"]))
        region_seq = ''.join([rsr.seq[abs_rel_map[s]] for s in range(31795202, 31795202+1)])
        print region_seq
        self.assertEqual(expected_seq, region_seq)


if __name__ == '__main__':
    unittest.main()
