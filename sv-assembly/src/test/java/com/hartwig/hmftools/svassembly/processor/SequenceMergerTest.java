package com.hartwig.hmftools.svassembly.processor;

import static org.assertj.core.api.Assertions.assertThat;

import com.hartwig.hmftools.svassembly.models.Sequence;

import org.junit.Test;

import htsjdk.samtools.SAMUtils;

public class SequenceMergerTest
{
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
        assertThat(merged.getBasesString()).isEqualTo(rb);
        assertThat(SAMUtils.phredToFastq(merged.getBaseQuality())).isEqualTo(rq);
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
        assertThat(merged.getBasesString()).isEqualTo("CGCGATCGTTGGTTGTCGAATCATATATATATATATATATTCGTCGTCG");
        assertThat(SAMUtils.phredToFastq(merged.getBaseQuality())).isEqualTo("F".repeat(lq.length()));
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
        assertThat(merged.getBasesString()).isEqualTo(expected);
        assertThat(SAMUtils.phredToFastq(merged.getBaseQuality())).isEqualTo("F".repeat(expected.length()));
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
        assertThat(merged.getBasesString()).isEqualTo(lb);
        assertThat(SAMUtils.phredToFastq(merged.getBaseQuality())).isEqualTo("F".repeat(lb.length()));
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
        assertThat(merged.getBasesString()).isEqualTo(expected);
        assertThat(SAMUtils.phredToFastq(merged.getBaseQuality())).isEqualTo("F".repeat(expected.length()));
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
        assertThat(merged.getBasesString()).isEqualTo(rb);
        assertThat(SAMUtils.phredToFastq(merged.getBaseQuality())).isEqualTo("F".repeat(rb.length()));
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
}
