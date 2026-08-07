package com.hartwig.hmftools.tars.common;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertFalse;
import static org.junit.Assert.assertNull;
import static org.junit.Assert.assertTrue;

import java.util.ArrayList;
import java.util.List;

import com.hartwig.hmftools.common.bam.CigarUtils;

import org.junit.Test;

import htsjdk.samtools.Cigar;
import htsjdk.samtools.CigarElement;
import htsjdk.samtools.CigarOperator;

public class TarsCigarUtilsTest
{
    @Test
    public void testNormalizeDropsZeroLengthAndMergesAdjacentOps()
    {
        assertEquals("72M2640N79M", TarsCigarUtils.normalize(cigar("72M102N0M2538N79M")).toString());
    }

    @Test
    public void testTerminalMatchedRunIgnoresTerminalSoftclipOnly()
    {
        List<CigarElement> elements = elements("5S30M10I4M6S");

        assertEquals(30, TarsCigarUtils.terminalMatchedRun(elements, false));
        assertEquals(4, TarsCigarUtils.terminalMatchedRun(elements, true));
        assertEquals(0, TarsCigarUtils.terminalMatchedRun(elements("5S2I30M"), false));
    }

    @Test
    public void testMatchPredicatesKeepMismatchOperatorExplicit()
    {
        assertTrue(TarsCigarUtils.isMatchedOp(CigarOperator.X));
        assertFalse(TarsCigarUtils.isMatchOrEqualOp(CigarOperator.X));
    }

    @Test
    public void testDetectsIndelNextToTerminalSoftclip()
    {
        assertTrue(TarsCigarUtils.indelAdjacentToTerminalSoftClip(elements("5S2I30M"), true));
        assertTrue(TarsCigarUtils.indelAdjacentToTerminalSoftClip(elements("30M2D5S"), false));
        assertFalse(TarsCigarUtils.indelAdjacentToTerminalSoftClip(elements("5S30M"), true));
        assertFalse(TarsCigarUtils.indelAdjacentToTerminalSoftClip(elements("30M5S"), false));
    }

    @Test
    public void testExtendTerminalSoftclipIntoAdjacentMatch()
    {
        List<CigarElement> leading = mutableElements("5S20M");
        assertTrue(TarsCigarUtils.extendTerminalSoftClipIntoMatch(leading, true, 3));
        assertEquals("2S23M", new Cigar(leading).toString());

        List<CigarElement> trailing = mutableElements("20M5S");
        assertTrue(TarsCigarUtils.extendTerminalSoftClipIntoMatch(trailing, false, 5));
        assertEquals("25M", new Cigar(trailing).toString());
    }

    @Test
    public void testRetractTerminalMatchIntoSoftclip()
    {
        assertEquals(
                "20M3M11S",
                new Cigar(TarsCigarUtils.retractTerminalMatchIntoSoftClip(elements("20M10M4S"), 7, true)).toString());
        assertEquals(
                "7S7M20M",
                new Cigar(TarsCigarUtils.retractTerminalMatchIntoSoftClip(elements("4S10M20M"), 3, false)).toString());
        assertNull(TarsCigarUtils.retractTerminalMatchIntoSoftClip(elements("4S3M20M"), 3, false));
    }

    @Test
    public void testClampReferenceOverhangIntoSoftclip()
    {
        assertEquals("1S9M", TarsCigarUtils.clampLeadingReferenceToSoftClip(cigar("10M"), 1).toString());
        assertEquals("9S78M", TarsCigarUtils.clampLeadingReferenceToSoftClip(cigar("6M1I80M"), 8).toString());
        assertEquals("40M4S", TarsCigarUtils.clampTrailingReferenceToSoftClip(cigar("40M1D4M"), 5).toString());
    }

    @Test
    public void testClampReturnsNullWhenOverhangCannotLeaveValidAlignedStart()
    {
        assertNull(TarsCigarUtils.clampLeadingReferenceToSoftClip(cigar("4M10S"), 6));
        assertNull(TarsCigarUtils.clampLeadingReferenceToSoftClip(cigar("3M1D80M"), 3));
    }

    private static Cigar cigar(final String cigar)
    {
        return CigarUtils.cigarFromStr(cigar);
    }

    private static List<CigarElement> elements(final String cigar)
    {
        return cigar(cigar).getCigarElements();
    }

    private static List<CigarElement> mutableElements(final String cigar)
    {
        return new ArrayList<>(elements(cigar));
    }
}
