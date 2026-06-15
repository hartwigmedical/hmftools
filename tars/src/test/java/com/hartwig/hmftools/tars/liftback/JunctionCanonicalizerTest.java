package com.hartwig.hmftools.tars.liftback;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertFalse;
import static org.junit.Assert.assertTrue;

import java.nio.charset.StandardCharsets;
import java.util.Arrays;

import org.junit.Test;

public class JunctionCanonicalizerTest
{
    private static final String CHR1 = "chr1";

    private static void put(final byte[] ref, final int pos, final String bases)
    {
        for(int i = 0; i < bases.length(); ++i)
            ref[pos - 1 + i] = (byte) bases.charAt(i);
    }

    private static JunctionCanonicalizer canonicalizer(final byte[] ref)
    {
        return new JunctionCanonicalizer(
                TarsTestFixtures.refSource(CHR1, new String(ref, StandardCharsets.US_ASCII)),
                JunctionCanonicalizer.DEFAULT_MAX_SHIFT);
    }

    // Shared fixture: read 5M10N5M starting at pos 1. The donor/acceptor at the bwa placement (shift 0)
    // are non-canonical (AA..CC); sliding the intron right by 2 lands it on GT..AG. The two read bases
    // that cross the intron under the +2 slide (read[5],read[6]) equal the donor-side ref (pos6,pos7).
    private static byte[] slidableRef()
    {
        final byte[] ref = new byte[25];
        Arrays.fill(ref, (byte) 'T');
        put(ref, 1, "CCCCC");  // left exon, matches read[0..4]
        put(ref, 6, "AA");     // donor at shift 0 -> non-canonical; also the +2 moved bases (match read)
        put(ref, 8, "GT");     // donor at shift +2 -> canonical
        put(ref, 14, "CC");    // acceptor at shift 0 -> non-canonical
        put(ref, 16, "AG");    // acceptor at shift +2 -> canonical
        return ref;
    }

    private static byte[] slidableRead()
    {
        // read[5]=read[6]='A' so the +2 slide's moved bases match ref pos6,pos7
        return "CCCCCAACCC".getBytes();
    }

    @Test
    public void slidesNonCanonicalJunctionToCanonical()
    {
        final JunctionCanonicalizationResult result = canonicalizer(slidableRef())
                .tryCanonicalize(CHR1, 1, "5M10N5M", slidableRead());

        assertTrue(result.Changed);
        assertEquals(1, result.JunctionsShifted);
        assertEquals("7M10N3M", result.NewCigar);
    }

    @Test
    public void leavesCanonicalJunctionUntouched()
    {
        // donor already GT at shift 0 (pos6), acceptor already AG (pos14): no slide attempted.
        final byte[] ref = new byte[25];
        Arrays.fill(ref, (byte) 'T');
        put(ref, 1, "CCCCC");
        put(ref, 6, "GT");
        put(ref, 14, "AG");

        final JunctionCanonicalizationResult result = canonicalizer(ref)
                .tryCanonicalize(CHR1, 1, "5M10N5M", "CCCCCTTCCC".getBytes());

        assertFalse(result.Changed);
    }

    @Test
    public void refusesSlideWhenMovedBasesMismatch()
    {
        // motif would become canonical at +2, but the read bases that must cross the intron (read[5,6])
        // do not match the donor-side ref (pos6,pos7 = AA), so the slide is unsafe and rejected.
        final byte[] ref = slidableRef();
        final byte[] read = "CCCCCGGCCC".getBytes(); // read[5,6] = GG != ref AA

        final JunctionCanonicalizationResult result = canonicalizer(ref)
                .tryCanonicalize(CHR1, 1, "5M10N5M", read);

        assertFalse(result.Changed);
    }

    @Test
    public void noRefSourceLeavesCigarUnchanged()
    {
        final JunctionCanonicalizer canon = new JunctionCanonicalizer(null, JunctionCanonicalizer.DEFAULT_MAX_SHIFT);
        assertFalse(canon.tryCanonicalize(CHR1, 1, "5M10N5M", slidableRead()).Changed);
    }
}
