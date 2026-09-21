package com.hartwig.hmftools.viridian.integration.seq_align;

import static com.hartwig.hmftools.common.genome.region.Orientation.FORWARD;

import static org.junit.Assert.assertEquals;

import com.hartwig.hmftools.viridian.reference.OncologyGroup;
import com.hartwig.hmftools.viridian.reference.ViralContig;

import org.junit.Test;

public class ViralSequenceAlignmentTest
{
    private static final ViralContig CONTIG = new ViralContig("v1", 7906, "Virus v1", new OncologyGroup("Group A"));

    // Clipped bases are the host side of an integration junction, so they are not insert bases placed on the virus.
    @Test
    public void excludesClippedBases()
    {
        assertEquals(70, alignment("30S70M", 100).alignedLength());
    }

    // Inserted query bases are not placed on the contig, and a deletion consumes no query bases at all.
    @Test
    public void excludesIndelBases()
    {
        assertEquals(28, alignment("10M2I8M3D10M", 30).alignedLength());
    }

    private static ViralSequenceAlignment alignment(String cigar, int sequenceLength)
    {
        return new ViralSequenceAlignment(CONTIG, 1500, FORWARD, cigar, 60, 1, sequenceLength);
    }
}
