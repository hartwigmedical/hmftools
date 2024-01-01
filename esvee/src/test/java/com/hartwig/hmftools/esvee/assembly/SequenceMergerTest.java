package com.hartwig.hmftools.esvee.assembly;

public class SequenceMergerTest
{
    /* CHASHA FIXME
    @Test
    public void canMergeCoincident()
    {
        final String lb = "ATCGTTGGTCGTCGAATC";
        final String lq = "FFFFFFFF!FFFFFFFFF";
        final String rb = "ATCGTTGGTTGTCGAATC";
        final String rq = "FFFFFFFFFFFFFFFFFF";
        final int supportIndex = 0;

        final Sequence left = Sequence.fromBytes(lb.getBytes(), SAMUtils.fastqToPhred(lq));
        final Sequence right = Sequence.fromBytes(rb.getBytes(), SAMUtils.fastqToPhred(rq));

        final Sequence merged = SequenceMerger.merge(left, right, supportIndex);
        assertTrue(merged.getBasesString().equals(rb));
        assertTrue(SAMUtils.phredToFastq(merged.getBaseQuality().equals(rq)));
    }

    @Test
    public void canMergeLengtheningRepeatLeft()
    {
        final String lb = "CGCGATCGTTGGTCGTCGAATCATATATATATATATATATTCGTCGTCG";
        final String lq = "FFFFFFFFFFFFF!FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF";
        final String rb =     "ATCGTTGGTTGTCGAATCATATATATATATATATATATTCGTCGTCG";
        final String rq =     "FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFCCFFFFFFFFFFF";
        final int supportIndex = 4;

        final Sequence left = Sequence.fromBytes(lb.getBytes(), SAMUtils.fastqToPhred(lq));
        final Sequence right = Sequence.fromBytes(rb.getBytes(), SAMUtils.fastqToPhred(rq));

        final Sequence merged = SequenceMerger.merge(left, right, supportIndex);
        assertTrue(merged.getBasesString().equals("CGCGATCGTTGGTTGTCGAATCATATATATATATATATATTCGTCGTCG"));
        assertTrue(SAMUtils.phredToFastq(merged.getBaseQuality()).equals("F".repeat(lq.length())));
    }

    @Test
    public void canMergeLengtheningRepeatRight()
    {
        final String lb = "CGCGATCGTTGGTCGTCGAATCATATATATATATATATATTCGTCGTCG";
        final String lq = "FFFFFFFFFFFFF!FFFFFFFFFFFFFFFFFFFFCCFFFFFFFFFFFFF";
        final String rb =     "ATCGTTGGTTGTCGAATCATATATATATATATATATATTCGTCGTCG";
        final String rq =     "FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF";
        final int supportIndex = 4;

        final Sequence left = Sequence.fromBytes(lb.getBytes(), SAMUtils.fastqToPhred(lq));
        final Sequence right = Sequence.fromBytes(rb.getBytes(), SAMUtils.fastqToPhred(rq));

        final Sequence merged = SequenceMerger.merge(left, right, supportIndex);
        final String expected = "CGCGATCGTTGGTTGTCGAATCATATATATATATATATATATTCGTCGTCG";
        assertTrue(merged.getBasesString().equals(expected));
        assertTrue(SAMUtils.phredToFastq(merged.getBaseQuality()).equals("F".repeat(expected.length())));
    }

    @Test
    public void canMergeLeftStart()
    {
        final String lb = "CGCGATCGTTGGTCGTCGAATCG";
        final String lq = "FFFFFFFF!FFFFFFFFFFFFFF";
        final String rb = "ATCGTTGGTTGTCGAATC";
        final String rq = "FFFFFFFFFEFFFFFFFF";
        final int supportIndex = 4;

        final Sequence left = Sequence.fromBytes(lb.getBytes(), SAMUtils.fastqToPhred(lq));
        final Sequence right = Sequence.fromBytes(rb.getBytes(), SAMUtils.fastqToPhred(rq));

        final Sequence merged = SequenceMerger.merge(left, right, supportIndex);
        assertTrue(merged.getBasesString().equals(lb));
        assertTrue(SAMUtils.phredToFastq(merged.getBaseQuality()).equals("F".repeat(lb.length())));
    }

    @Test
    public void canMergeLeftoverLeft()
    {
        final String lb = "ATCGTTGGTCGTCGAATCGCGC";
        final String lq = "FFFFFFFF!FFFFFFFFFFFFF";
        final String rb = "ATCGTTGGTTGTCGAATC";
        final String rq = "FFFFFFFFFFFFFFFFFF";
        final int supportIndex = 0;

        final Sequence left = Sequence.fromBytes(lb.getBytes(), SAMUtils.fastqToPhred(lq));
        final Sequence right = Sequence.fromBytes(rb.getBytes(), SAMUtils.fastqToPhred(rq));

        final Sequence merged = SequenceMerger.merge(left, right, supportIndex);
        final String expected = "ATCGTTGGTTGTCGAATCGCGC";
        assertTrue(merged.getBasesString().equals(expected));
        assertTrue(SAMUtils.phredToFastq(merged.getBaseQuality()).equals("F".repeat(expected.length())));
    }

    @Test
    public void canMergeLeftoverRight()
    {
        final String lb = "ATCGTTGGTCGTCGAATC";
        final String lq = "FFFFFFFFF!FFFFFFFF";
        final String rb = "ATCGTTGGTTGTCGAATCCGCG";
        final String rq = "FFFFFFFFFFFFFFFFFFFFFF";
        final int supportIndex = 0;

        final Sequence left = Sequence.fromBytes(lb.getBytes(), SAMUtils.fastqToPhred(lq));
        final Sequence right = Sequence.fromBytes(rb.getBytes(), SAMUtils.fastqToPhred(rq));

        final Sequence merged = SequenceMerger.merge(left, right, supportIndex);
        assertTrue(merged.getBasesString().equals(rb));
        assertTrue(SAMUtils.phredToFastq(merged.getBaseQuality()).equals("F".repeat(rb.length())));
    }

    @Test
    public void canMergeComplex()
    {
        final String lb = "AAGTAATTGCGGTTTTTGCCAAATGGCAAAAACCGCAATTACTTGTGCACCAACCTAATACTATAATTCACCAGATGCATTATAGTATGAATGAAGAGAC";
        final String lq = "F".repeat(lb.length());
        final String rb = "AAGTAATTTCGGTTTTTTGCAAATGGCAAAAACCGCAATTACTTGTGCACCAACCTAATACTATAATTCACCAGATGCATTATAGTATGAATGAAGAGAC";
        final String rq = "F".repeat(rb.length());
        final int supportIndex = 0;

        final Sequence left = Sequence.fromBytes(lb.getBytes(), SAMUtils.fastqToPhred(lq));
        final Sequence right = Sequence.fromBytes(rb.getBytes(), SAMUtils.fastqToPhred(rq));

        final Sequence merged = SequenceMerger.merge(left, right, supportIndex);

    }
     */
}
