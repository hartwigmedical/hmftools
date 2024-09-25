package com.hartwig.hmftools.sage.common;

import static com.hartwig.hmftools.sage.common.UltimaVariantReadContextBuilderUtils.isAdjacentToLongHomopolymer;
import static com.hartwig.hmftools.sage.common.UltimaVariantReadContextBuilderUtils.isMsiIndelOfType;

import static org.junit.Assert.assertFalse;
import static org.junit.Assert.assertTrue;
import static org.mockito.Mockito.mock;
import static org.mockito.Mockito.when;

import java.util.List;

import com.google.common.collect.Lists;

import org.junit.Test;

import htsjdk.samtools.SAMRecord;

public class UltimaVariantReadContextBuilderUtilsTest
{
    @Test
    public void testIsMsiIndelOfTypeNotIndel()
    {
        SimpleVariant mockVariant = mock(SimpleVariant.class);
        when(mockVariant.isIndel()).thenReturn(false);

        assertFalse(isMsiIndelOfType(mockVariant, null));
    }

    @Test
    public void testIsMsiIndelOfTypeNoUnits()
    {
        SimpleVariant mockVariant = mock(SimpleVariant.class);
        when(mockVariant.isIndel()).thenReturn(true);
        when(mockVariant.isInsert()).thenReturn(true);
        when(mockVariant.ref()).thenReturn("A");
        when(mockVariant.alt()).thenReturn("AA");

        List<String> units = Lists.newArrayList();
        assertFalse(isMsiIndelOfType(mockVariant, units));
    }

    @Test
    public void testIsMsiIndelOfTypeSingleUnitNotMsiIndel()
    {
        SimpleVariant mockVariant = mock(SimpleVariant.class);
        when(mockVariant.isIndel()).thenReturn(true);
        when(mockVariant.isInsert()).thenReturn(true);
        when(mockVariant.ref()).thenReturn("A");
        when(mockVariant.alt()).thenReturn("AAB");

        List<String> units = Lists.newArrayList("A");
        assertFalse(isMsiIndelOfType(mockVariant, units));
    }

    @Test
    public void testIsMsiIndelOfTypeSingleUnitIsMsiIndel()
    {
        SimpleVariant mockVariant = mock(SimpleVariant.class);
        when(mockVariant.isIndel()).thenReturn(true);
        when(mockVariant.isInsert()).thenReturn(true);
        when(mockVariant.ref()).thenReturn("A");
        when(mockVariant.alt()).thenReturn("AABABAB");

        List<String> units = Lists.newArrayList("AB");
        assertTrue(isMsiIndelOfType(mockVariant, units));
    }

    @Test
    public void testIsAdjacentToLongHomopolymerNotLongEnough()
    {
        int longLength = 10;
        SAMRecord mockRead = mock(SAMRecord.class);
        when(mockRead.getReadBases()).thenReturn("A".repeat(2*longLength - 1).getBytes());

        assertFalse(isAdjacentToLongHomopolymer(mockRead, longLength - 1, longLength));
    }

    @Test
    public void testIsAdjacentToLongHomopolymerOnLeft()
    {
        int longLength = 10;
        String readBases = "A".repeat(longLength) + "C" + "T".repeat(longLength - 1);
        SAMRecord mockRead = mock(SAMRecord.class);
        when(mockRead.getReadBases()).thenReturn(readBases.getBytes());

        assertTrue(isAdjacentToLongHomopolymer(mockRead, longLength, longLength));
    }

    @Test
    public void testIsAdjacentToLongHomopolymerOnRight()
    {
        int longLength = 10;
        String readBases = "T".repeat(longLength - 1) + "C" + "A".repeat(longLength);
        SAMRecord mockRead = mock(SAMRecord.class);
        when(mockRead.getReadBases()).thenReturn(readBases.getBytes());

        assertTrue(isAdjacentToLongHomopolymer(mockRead, longLength - 1, longLength));
    }

    @Test
    public void testIsAdjacentToLongHomopolymerNotOnLeft()
    {
        int longLength = 10;
        String readBases = "C" + "A".repeat(2*longLength - 1);
        SAMRecord mockRead = mock(SAMRecord.class);
        when(mockRead.getReadBases()).thenReturn(readBases.getBytes());

        assertFalse(isAdjacentToLongHomopolymer(mockRead, longLength, longLength));
    }

    @Test
    public void testIsAdjacentToLongHomopolymerNotOnRight()
    {
        int longLength = 10;
        String readBases = "A".repeat(2*longLength - 1) + "C";
        SAMRecord mockRead = mock(SAMRecord.class);
        when(mockRead.getReadBases()).thenReturn(readBases.getBytes());

        assertFalse(isAdjacentToLongHomopolymer(mockRead, longLength - 1, longLength));
    }

    @Test
    public void testIsAdjacentToLongHomopolymerNHomopolymerOnLeft()
    {
        int longLength = 10;
        SAMRecord mockRead = mock(SAMRecord.class);
        when(mockRead.getReadBases()).thenReturn("N".repeat(2*longLength).getBytes());

        assertFalse(isAdjacentToLongHomopolymer(mockRead, longLength, longLength));
    }
}
