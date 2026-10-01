package com.hartwig.hmftools.viridian.integration.align;

import static com.hartwig.hmftools.common.bam.CigarUtils.cigarFromStr;
import static com.hartwig.hmftools.common.genome.region.Orientation.FORWARD;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertFalse;
import static org.junit.Assert.assertTrue;

import com.hartwig.hmftools.viridian.reference.OncologyGroup;
import com.hartwig.hmftools.viridian.reference.ViralContig;

import org.junit.Test;

public class ViralInsertAlignmentTest
{
    private static final ViralContig CONTIG = new ViralContig("v1", 7906, "Virus v1", new OncologyGroup("Group A"));

    // Clipped bases are the host side of an integration junction, so they are not insert bases placed on the virus.
    @Test
    public void testAlignedLengthExcludesClippedBases()
    {
        assertEquals(70, alignment("30S70M", 100).alignedLength());
    }

    // Inserted query bases are not placed on the contig, and a deletion consumes no query bases at all.
    @Test
    public void testAlignedLengthExcludesIndelBases()
    {
        assertEquals(28, alignment("10M2I8M3D10M", 30).alignedLength());
    }

    @Test
    public void testScorePerAlignedBase()
    {
        assertEquals(60.0 / 70, alignment("30S70M", 100, 60).scorePerAlignedBase(), 1e-9);
    }

    @Test
    public void testPassesFiltersWhenScoreAndPerBaseClearFloors()
    {
        assertTrue(alignment("70M", 70, 60).passesFilters());
    }

    @Test
    public void testDoesNotClearWhenScoreBelowFloor()
    {
        assertFalse(alignment("70M", 70, 25).passesFilters());
    }

    // A long match reaching the score floor on length alone, at low identity, is rejected by the per-aligned-base floor.
    @Test
    public void testDoesNotClearWhenPerBaseBelowFloor()
    {
        assertFalse(alignment("100M", 100, 40).passesFilters());
    }

    @Test
    public void testClearsAtPerBaseFloor()
    {
        assertTrue(alignment("50M", 50, 35).passesFilters());
    }

    private static ViralInsertAlignment alignment(String cigar, int sequenceLength)
    {
        return alignment(cigar, sequenceLength, 60);
    }

    private static ViralInsertAlignment alignment(String cigar, int sequenceLength, int alignerScore)
    {
        return new ViralInsertAlignment(CONTIG, 1500, FORWARD, cigarFromStr(cigar), alignerScore, 1, sequenceLength);
    }
}
