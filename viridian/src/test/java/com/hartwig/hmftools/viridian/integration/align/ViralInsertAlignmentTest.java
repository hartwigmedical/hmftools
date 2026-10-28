package com.hartwig.hmftools.viridian.integration.align;

import static com.hartwig.hmftools.common.bam.CigarUtils.cigarFromStr;
import static com.hartwig.hmftools.common.genome.region.Orientation.FORWARD;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertFalse;
import static org.junit.Assert.assertTrue;

import com.hartwig.hmftools.viridian.reference.OncologyGroup;
import com.hartwig.hmftools.viridian.reference.VirusGenome;

import org.junit.Test;

public class ViralInsertAlignmentTest
{
    private static final VirusGenome GENOME = new VirusGenome("v1", 7906, "Virus v1", new OncologyGroup("Group A"));

    // Clipped bases are the host side of an integration junction, so they are not insert bases placed on the virus.
    @Test
    public void testAlignedLengthExcludesClippedBases()
    {
        assertEquals(70, alignment("30S70M", 100).alignedLength());
    }

    @Test
    public void testAlignedLengthOfPlainMatch()
    {
        assertEquals(70, alignment("70M", 70).alignedLength());
    }

    @Test
    public void testAlignedLengthCountsInternalInsertion()
    {
        assertEquals(20, alignment("10M2I8M", 20).alignedLength());
    }

    @Test
    public void testAlignedLengthExcludesDeletion()
    {
        assertEquals(20, alignment("10M3D10M", 20).alignedLength());
    }

    @Test
    public void testAlignedLengthExcludesHardClippedBases()
    {
        assertEquals(70, alignment("30H70M", 100).alignedLength());
    }

    @Test
    public void testPassesFiltersWhenScoreClearsFloor()
    {
        assertTrue(alignment("70M", 70, 60).passesFilters());
    }

    @Test
    public void testDoesNotClearWhenScoreBelowFloor()
    {
        assertFalse(alignment("70M", 70, 25).passesFilters());
    }

    // A long, low-identity match (score 40 over 100 bases) passes on score alone.
    @Test
    public void testPassesFiltersAtLowScorePerAlignedBase()
    {
        assertTrue(alignment("100M", 100, 40).passesFilters());
    }

    private static ViralInsertAlignment alignment(String cigar, int sequenceLength)
    {
        return alignment(cigar, sequenceLength, 60);
    }

    private static ViralInsertAlignment alignment(String cigar, int sequenceLength, int alignerScore)
    {
        return new ViralInsertAlignment(GENOME, 1500, FORWARD, cigarFromStr(cigar), alignerScore, 1, sequenceLength);
    }
}
