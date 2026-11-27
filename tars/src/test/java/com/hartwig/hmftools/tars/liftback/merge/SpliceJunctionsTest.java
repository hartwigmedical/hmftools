package com.hartwig.hmftools.tars.liftback.merge;

import static com.hartwig.hmftools.common.test.GeneTestUtils.CHR_1;
import static com.hartwig.hmftools.tars.liftback.TarsTestFixtures.bases;
import static com.hartwig.hmftools.tars.liftback.TarsTestFixtures.refGenome;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertNull;

import java.util.List;

import com.hartwig.hmftools.common.genome.region.Orientation;
import com.hartwig.hmftools.common.region.BaseRegion;
import com.hartwig.hmftools.tars.common.ContigEntry;
import com.hartwig.hmftools.tars.liftback.EnsemblAnnotationIndex;
import com.hartwig.hmftools.tars.liftback.merge.SpliceJunctions.Tier;

import org.junit.Test;

public class SpliceJunctionsTest
{
    @Test
    public void testClassifiesMotifTiers()
    {
        // Strand is unknown at scan time, so GT-AG and its reverse complement CT-AC both score canonical.
        assertEquals(Tier.CANONICAL, SpliceJunctions.motifTier(bases("GT"), bases("AG")));
        assertEquals(Tier.CANONICAL, SpliceJunctions.motifTier(bases("CT"), bases("AC")));

        assertEquals(Tier.SEMI_CANONICAL, SpliceJunctions.motifTier(bases("GC"), bases("AG")));
        assertEquals(Tier.SEMI_CANONICAL, SpliceJunctions.motifTier(bases("CT"), bases("GC")));
        assertEquals(Tier.SEMI_CANONICAL, SpliceJunctions.motifTier(bases("AT"), bases("AC")));
        assertEquals(Tier.SEMI_CANONICAL, SpliceJunctions.motifTier(bases("GT"), bases("AT")));

        assertEquals(Tier.CANONICAL, SpliceJunctions.motifTier(bases("gt"), bases("ag")));
        assertEquals(Tier.CANONICAL, SpliceJunctions.motifTier(bases("ct"), bases("ac")));

        // Both flanks must match: a canonical donor with a non-motif acceptor is NONE.
        assertEquals(Tier.NONE, SpliceJunctions.motifTier(bases("AA"), bases("GG")));
        assertEquals(Tier.NONE, SpliceJunctions.motifTier(bases("NN"), bases("NN")));
        assertEquals(Tier.NONE, SpliceJunctions.motifTier(bases("GT"), bases("CC")));
        assertEquals(Tier.NONE, SpliceJunctions.motifTier(bases("CC"), bases("AG")));
    }

    @Test
    public void testMotifStrandFollowsTheOrientedMotif()
    {
        assertEquals(Orientation.FORWARD, SpliceJunctions.motifStrand(bases("GT"), bases("AG")));
        assertEquals(Orientation.REVERSE, SpliceJunctions.motifStrand(bases("CT"), bases("AC")));
        assertNull(SpliceJunctions.motifStrand(bases("AA"), bases("GG")));
    }

    @Test
    public void testSpliceStrandComesFromTheAnnotatedTranscriptStrand()
    {
        SpliceJunctions spliceJunctions = new SpliceJunctions(reverseStrandTranscript(), refGenome(CHR_1, "A".repeat(500)));

        assertEquals(Orientation.REVERSE, spliceJunctions.spliceStrand(CHR_1, 150, "50M100N50M"));

        // no N gap, so no junction to take a strand from
        assertNull(spliceJunctions.spliceStrand(CHR_1, 150, "100M"));

        // an N gap that is not an annotated junction, with no ref bases to fall back on
        assertNull(spliceJunctions.spliceStrand(CHR_1, 150, "50M150N50M"));
    }

    private static EnsemblAnnotationIndex reverseStrandTranscript()
    {
        return EnsemblAnnotationIndex.fromContigEntries(List.of(
                ContigEntry.annotationOnly(
                        "g", "gn", "tn", CHR_1, -1,
                        List.of(new BaseRegion(100, 199), new BaseRegion(300, 399)))));
    }
}
