package com.hartwig.hmftools.tars.liftback;

import static com.hartwig.hmftools.tars.liftback.TarsTestFixtures.bases;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertFalse;
import static org.junit.Assert.assertTrue;

import com.hartwig.hmftools.tars.common.TarsConstants;
import com.hartwig.hmftools.tars.liftback.TarsTestFixtures.TestGenome;

import org.junit.Test;

public class JunctionCanonicalizerTest
{
    private static final String CHR1 = "chr1";

    private static JunctionCanonicalizer canonicalizer(final TestGenome genome)
    {
        return new JunctionCanonicalizer(genome.asRefSource(), TarsConstants.DEFAULT_MAX_SHIFT);
    }

    // Shared fixture: read 5M10N5M starting at pos 1. The donor/acceptor at the bwa placement (shift 0)
    // are non-canonical (AA..CC); sliding the intron right by 2 lands it on GT..AG. The two read bases
    // that cross the intron under the +2 slide (read[5],read[6]) equal the donor-side ref (pos6,pos7).
    private static TestGenome slidableGenome()
    {
        return new TestGenome().with(CHR1, 25, 'T')
                .set(CHR1, 1, "CCCCC")  // left exon, matches read[0..4]
                .set(CHR1, 6, "AA")     // donor at shift 0 -> non-canonical; also the +2 moved bases (match read)
                .set(CHR1, 8, "GT")     // donor at shift +2 -> canonical
                .set(CHR1, 14, "CC")    // acceptor at shift 0 -> non-canonical
                .set(CHR1, 16, "AG");   // acceptor at shift +2 -> canonical
    }

    private static byte[] slidableRead()
    {
        // read[5]=read[6]='A' so the +2 slide's moved bases match ref pos6,pos7
        return bases("CCCCCAACCC");
    }

    @Test
    public void slidesNonCanonicalJunctionToCanonical()
    {
        JunctionCanonicalizer canonicalizer = canonicalizer(slidableGenome());
        JunctionCanonicalizationResult result = canonicalizer.tryCanonicalize(CHR1, 1, "5M10N5M", slidableRead());

        assertTrue(result.changed());
        assertEquals(1, canonicalizer.junctionsShifted());
        assertEquals("7M10N3M", result.newCigar());
    }

    @Test
    public void leavesCanonicalJunctionUntouched()
    {
        // donor already GT at shift 0 (pos6), acceptor already AG (pos14): no slide attempted.
        TestGenome genome = new TestGenome().with(CHR1, 25, 'T')
                .set(CHR1, 1, "CCCCC").set(CHR1, 6, "GT").set(CHR1, 14, "AG");

        JunctionCanonicalizationResult result = canonicalizer(genome)
                .tryCanonicalize(CHR1, 1, "5M10N5M", bases("CCCCCTTCCC"));

        assertFalse(result.changed());
    }

    @Test
    public void refusesSlideWhenMovedBasesMismatch()
    {
        // motif would become canonical at +2, but the read bases that must cross the intron (read[5,6])
        // do not match the donor-side ref (pos6,pos7 = AA), so the slide is unsafe and rejected.
        JunctionCanonicalizationResult result = canonicalizer(slidableGenome())
                .tryCanonicalize(CHR1, 1, "5M10N5M", bases("CCCCCGGCCC"));   // read[5,6] = GG != ref AA

        assertFalse(result.changed());
    }

    @Test
    public void noRefSourceLeavesCigarUnchanged()
    {
        JunctionCanonicalizer canon = new JunctionCanonicalizer(null, TarsConstants.DEFAULT_MAX_SHIFT);
        assertFalse(canon.tryCanonicalize(CHR1, 1, "5M10N5M", slidableRead()).changed());
    }
}
