package com.hartwig.hmftools.esvee.common.saga;

import static com.hartwig.hmftools.common.bam.CigarUtils.cigarFromStr;
import static com.hartwig.hmftools.common.genome.region.Orientation.FORWARD;
import static com.hartwig.hmftools.common.genome.region.Orientation.REVERSE;

import static org.junit.Assert.assertEquals;

import java.util.List;

import com.hartwig.hmftools.common.region.BasePosition;

import org.broadinstitute.hellbender.utils.bwa.BwaMemAlignment;
import org.junit.Test;

import htsjdk.samtools.Cigar;

public class SagaAlignmentTest
{
    private static SagaAssembly assembly(int length)
    {
        // Minimal valid deletion assembly with specified sequence length.
        String sequence = "A".repeat(length);
        SagaVariant variant = new SagaVariant(
                "Test",
                new SagaBreakend(new BasePosition("chr1", 100), FORWARD),
                new SagaBreakend(new BasePosition("chr1", 200), REVERSE),
                ""
        );
        return new SagaAssembly("Test", variant, List.of(length / 2), sequence);
    }

    private static BwaMemAlignment rawAlignment(int refStart, int refEnd)
    {
        return new BwaMemAlignment(0, 0, refStart, refEnd, 0, 0, 60, 0, 0, 0, "10M", "", "", 0, 0, 0);
    }

    private static SagaAlignment alignment(int queryLength, int sagaStart, int sagaEnd, final String cigarStr, int sagaSequenceLength)
    {
        Cigar cigar = cigarFromStr(cigarStr);
        BwaMemAlignment raw = rawAlignment(sagaStart, sagaEnd);
        SagaAssembly sagaAssembly = assembly(sagaSequenceLength);
        return new SagaAlignment(raw, cigar, queryLength, sagaAssembly);
    }

    @Test
    public void testLeftUnalignedQueryStartLimits()
    {
        // queryStart=5 (5S clip), sagaStart=10 → min(5, 10) = 5
        SagaAlignment a = alignment(25, 10, 25, "5S15M5S", 60);
        assertEquals(5, a.leftUnaligned());
    }

    @Test
    public void testLeftUnalignedSagaStartLimits()
    {
        // queryStart=10 (10S clip), sagaStart=3 → min(10, 3) = 3
        SagaAlignment a = alignment(25, 3, 15, "10S12M3S", 60);
        assertEquals(3, a.leftUnaligned());
    }

    @Test
    public void testLeftUnalignedBothEqual()
    {
        // queryStart=5, sagaStart=5 → min(5, 5) = 5
        SagaAlignment a = alignment(20, 5, 15, "5S10M5S", 60);
        assertEquals(5, a.leftUnaligned());
    }

    @Test
    public void testLeftUnalignedNoClipNoSagaOffset()
    {
        // queryStart=0 (no clip), sagaStart=0 → min(0, 0) = 0
        SagaAlignment a = alignment(10, 0, 10, "10M", 20);
        assertEquals(0, a.leftUnaligned());
    }

    @Test
    public void testRightUnalignedQuerySideLimits()
    {
        // rightClip=5, sagaLength - sagaEnd = 60-25=35 → min(5, 35) = 5
        SagaAlignment a = alignment(25, 10, 20, "5S10M5S", 60);
        assertEquals(5, a.rightUnaligned());
    }

    @Test
    public void testRightUnalignedSagaSideLimits()
    {
        // rightClip=10, sagaLength - sagaEnd = 20-15=5 → min(10, 5) = 5
        SagaAlignment a = alignment(25, 5, 15, "5S10M10S", 20);
        assertEquals(5, a.rightUnaligned());
    }

    @Test
    public void testRightUnalignedBothEqual()
    {
        // rightClip=5, sagaLength - sagaEnd = 25-20=5 → min(5, 5) = 5
        SagaAlignment a = alignment(20, 10, 20, "5S10M5S", 25);
        assertEquals(5, a.rightUnaligned());
    }

    @Test
    public void testRightUnalignedNoClipNoSagaRemaining()
    {
        // rightClip=0, sagaLength - sagaEnd = 10-10=0 → min(0, 0) = 0
        SagaAlignment a = alignment(10, 0, 10, "10M", 10);
        assertEquals(0, a.rightUnaligned());
    }

    @Test
    public void testSagaIndexToQueryIndexSimpleMatch()
    {
        // Cigar "10M", sagaStart=5, queryStart=0 (no clip)
        // sagaIndex=5 → queryIndex=0, sagaIndex=10 → queryIndex=5, sagaIndex=14 → queryIndex=9
        SagaAlignment a = alignment(10, 5, 15, "10M", 30);
        assertEquals(0, a.sagaIndexToQueryIndex(5));
        assertEquals(5, a.sagaIndexToQueryIndex(10));
        assertEquals(9, a.sagaIndexToQueryIndex(14));
    }

    @Test
    public void testSagaIndexToQueryIndexLeftExtrapolationBelowSagaStart()
    {
        // Cigar "5S10M", sagaStart=10, queryStart=5
        // sagaIndex=8 (< sagaStart): queryStart - (sagaStart - sagaIndex) = 5 - 2 = 3
        SagaAlignment a = alignment(15, 10, 20, "5S10M", 40);
        assertEquals(3, a.sagaIndexToQueryIndex(8));
    }

    @Test
    public void testSagaIndexToQueryIndexLeftExtrapolationAtSagaStart()
    {
        // sagaIndex=sagaStart: extrapolation returns queryStart
        SagaAlignment a = alignment(15, 10, 20, "5S10M", 40);
        assertEquals(5, a.sagaIndexToQueryIndex(10));
    }

    @Test
    public void testSagaIndexToQueryIndexRightExtrapolationAboveSagaEnd()
    {
        // Cigar "10M5S", sagaStart=10, sagaEnd=20, queryEnd=10 (queryLength=15, rightClip=5)
        // sagaIndex=22 (> sagaEnd): queryEnd + (22 - sagaEnd) = 10 + 2 = 12
        SagaAlignment a = alignment(15, 10, 20, "10M5S", 40);
        assertEquals(12, a.sagaIndexToQueryIndex(22));
    }

    @Test
    public void testSagaIndexToQueryIndexRightExtrapolationAtSagaEnd()
    {
        // sagaIndex=sagaEnd: extrapolation returns queryEnd
        SagaAlignment a = alignment(15, 10, 20, "10M5S", 40);
        assertEquals(10, a.sagaIndexToQueryIndex(20));
    }

    @Test
    public void testSagaIndexToQueryIndexDeletion()
    {
        // Cigar "5M3D5M": sagaStart=10, sagaEnd=23, queryLength=10
        // Positions in D (sagaStart+5=15, 16, 17) map to next read base = index 5
        SagaAlignment a = alignment(10, 10, 23, "5M3D5M", 40);
        assertEquals(5, a.sagaIndexToQueryIndex(15));
        assertEquals(5, a.sagaIndexToQueryIndex(17));
    }

    @Test
    public void testSagaIndexToQueryIndexInsertion()
    {
        // Cigar "5M3I5M": sagaStart=10, sagaEnd=20, queryLength=13
        // sagaIndex=13 (within first 5M): index = 3
        // sagaIndex=15 (= sagaStart+5, start of second M block after 3I): index = 8 (5 ref-bases + 3 inserted)
        SagaAlignment a = alignment(13, 10, 20, "5M3I5M", 40);
        assertEquals(3, a.sagaIndexToQueryIndex(13));
        assertEquals(8, a.sagaIndexToQueryIndex(15));
    }
}
