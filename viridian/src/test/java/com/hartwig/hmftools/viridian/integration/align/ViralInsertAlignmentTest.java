package com.hartwig.hmftools.viridian.integration.align;

import static com.hartwig.hmftools.common.bam.CigarUtils.cigarFromStr;
import static com.hartwig.hmftools.common.genome.region.Orientation.FORWARD;
import static com.hartwig.hmftools.viridian.common.ViridianConstants.INTEGRATION_ALIGN_SCORE_MIN;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertFalse;
import static org.junit.Assert.assertTrue;

import com.hartwig.hmftools.viridian.reference.Virus;
import com.hartwig.hmftools.viridian.reference.VirusGenome;
import com.hartwig.hmftools.viridian.reference.VirusGroup;

import org.junit.Test;

public class ViralInsertAlignmentTest
{
    private static final VirusGenome GENOME =
            new VirusGenome("v1", 7906, new Virus("Virus v1", new VirusGroup("Group A")));

    @Test
    public void testAlignedLength()
    {
        // Counts the insert bases placed on the virus: matched and inserted bases. Clipped bases are the host side of an
        // integration junction, and deleted bases are virus bases, not insert bases.

        assertEquals(30, alignment("10S10M2I8M3D10M7S", 47).alignedLength());
        assertEquals(30, alignment("10H10M2I8M3D10M7H", 47).alignedLength());
    }

    @Test
    public void testPassesFiltersWhenScoreClearsFloor()
    {
        assertTrue(alignment("70M", 70, INTEGRATION_ALIGN_SCORE_MIN).passesFilters());
    }

    @Test
    public void testPassesFiltersFailsWhenScoreBelowFloor()
    {
        assertFalse(alignment("70M", 70, INTEGRATION_ALIGN_SCORE_MIN - 1).passesFilters());
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
