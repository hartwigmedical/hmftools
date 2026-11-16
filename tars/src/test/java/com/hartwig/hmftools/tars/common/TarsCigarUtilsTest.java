package com.hartwig.hmftools.tars.common;

import static com.hartwig.hmftools.common.bam.CigarUtils.cigarElementsFromStr;
import static com.hartwig.hmftools.common.bam.CigarUtils.cigarFromStr;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertFalse;
import static org.junit.Assert.assertNull;
import static org.junit.Assert.assertTrue;

import java.util.List;

import org.junit.Test;

import htsjdk.samtools.CigarElement;
import htsjdk.samtools.CigarOperator;

public class TarsCigarUtilsTest
{
    @Test
    public void testNormalizeDropsZeroLengthAndMergesAdjacentOps()
    {
        assertEquals("72M2640N79M", TarsCigarUtils.normalize(cigarFromStr("72M102N0M2538N79M")).toString());
    }

    @Test
    public void testTerminalMatchedRunIgnoresTerminalSoftclipOnly()
    {
        List<CigarElement> elements = cigarElementsFromStr("5S30M10I4M6S");

        assertEquals(30, TarsCigarUtils.terminalMatchedRun(elements, false));
        assertEquals(4, TarsCigarUtils.terminalMatchedRun(elements, true));
        assertEquals(0, TarsCigarUtils.terminalMatchedRun(cigarElementsFromStr("5S2I30M"), false));
        assertEquals(30, TarsCigarUtils.terminalMatchedRun(cigarElementsFromStr("5S30X10I4M"), false));
    }

    @Test
    public void testMatchOrEqualExcludesMismatchOperator()
    {
        assertFalse(TarsCigarUtils.isMatchOrEqualOp(CigarOperator.X));
    }

    @Test
    public void testDetectsIndelNextToTerminalSoftclip()
    {
        assertTrue(TarsCigarUtils.indelAdjacentToTerminalSoftClip(cigarElementsFromStr("5S2I30M"), true));
        assertTrue(TarsCigarUtils.indelAdjacentToTerminalSoftClip(cigarElementsFromStr("30M2D5S"), false));
        assertFalse(TarsCigarUtils.indelAdjacentToTerminalSoftClip(cigarElementsFromStr("5S30M"), true));
        assertFalse(TarsCigarUtils.indelAdjacentToTerminalSoftClip(cigarElementsFromStr("30M5S"), false));
    }

    @Test
    public void testClampReferenceOverhangIntoSoftclip()
    {
        assertEquals("1S9M", TarsCigarUtils.clampLeadingReferenceToSoftClip(cigarFromStr("10M"), 1).toString());
        assertEquals("9S78M", TarsCigarUtils.clampLeadingReferenceToSoftClip(cigarFromStr("6M1I80M"), 8).toString());
        assertEquals("40M4S", TarsCigarUtils.clampTrailingReferenceToSoftClip(cigarFromStr("40M1D4M"), 5).toString());
    }

    @Test
    public void testClampReturnsNullWhenOverhangCannotLeaveValidAlignedStart()
    {
        assertNull(TarsCigarUtils.clampLeadingReferenceToSoftClip(cigarFromStr("4M10S"), 6));
        assertNull(TarsCigarUtils.clampLeadingReferenceToSoftClip(cigarFromStr("3M1D80M"), 3));
    }

}
