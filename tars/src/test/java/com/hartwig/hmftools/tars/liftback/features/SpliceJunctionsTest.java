package com.hartwig.hmftools.tars.liftback.features;

import static com.hartwig.hmftools.common.test.GeneTestUtils.CHR_1;

import static org.junit.Assert.assertEquals;

import java.nio.charset.StandardCharsets;
import java.util.List;

import com.hartwig.hmftools.common.region.BaseRegion;
import com.hartwig.hmftools.tars.common.ContigEntry;
import com.hartwig.hmftools.tars.liftback.EnsemblAnnotationIndex;
import com.hartwig.hmftools.tars.liftback.features.SpliceJunctions.Tier;

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
        assertEquals(1, SpliceJunctions.motifStrand(bases("GT"), bases("AG")));
        assertEquals(-1, SpliceJunctions.motifStrand(bases("CT"), bases("AC")));
        assertEquals(0, SpliceJunctions.motifStrand(bases("AA"), bases("GG")));
    }

    @Test
    public void testRejectsNullOrWrongLengthInputs()
    {
        assertEquals(Tier.NONE, SpliceJunctions.motifTier(null, bases("AG")));
        assertEquals(Tier.NONE, SpliceJunctions.motifTier(bases("GT"), null));
        assertEquals(Tier.NONE, SpliceJunctions.motifTier(bases("G"), bases("AG")));
        assertEquals(Tier.NONE, SpliceJunctions.motifTier(bases("GT"), bases("A")));
        assertEquals(Tier.NONE, SpliceJunctions.motifTier(bases("GTC"), bases("AG")));
    }

    @Test
    public void testSpliceStrandComesFromTheAnnotatedTranscriptStrand()
    {
        SpliceJunctions spliceJunctions = new SpliceJunctions(reverseStrandTranscript(), null);

        assertEquals(-1, spliceJunctions.spliceStrand(CHR_1, 150, "50M100N50M"));

        // no N gap, so no junction to take a strand from
        assertEquals(0, spliceJunctions.spliceStrand(CHR_1, 150, "100M"));

        // an N gap that is not an annotated junction, with no ref genome to fall back on
        assertEquals(0, spliceJunctions.spliceStrand(CHR_1, 150, "50M150N50M"));
    }

    private static EnsemblAnnotationIndex reverseStrandTranscript()
    {
        return EnsemblAnnotationIndex.fromContigEntries(List.of(
                ContigEntry.annotationOnly(
                        "g", "gn", "tn", CHR_1, -1,
                        List.of(new BaseRegion(100, 199), new BaseRegion(300, 399)))));
    }

    private static byte[] bases(final String bases)
    {
        return bases.getBytes(StandardCharsets.US_ASCII);
    }
}
