package com.hartwig.hmftools.viridian.common;

import static com.hartwig.hmftools.viridian.common.Utils.parseMdTagMismatchPositions;

import static org.junit.Assert.assertArrayEquals;
import static org.junit.Assert.assertThrows;

import org.junit.Test;

public class UtilsTest
{
    @Test
    public void testParseMdTagMismatchPositionsSkipsMatchesAndDeletions()
    {
        // Mismatch at 100, matches 101-102, mismatch at 103, matches 104-115, deletion 116-117,
        // mismatch at 118 straight after the deletion, matches 119-122.

        assertArrayEquals(new int[] { 100, 103, 118 }, parseMdTagMismatchPositions(100, "0A2T12^GC0C4"));
    }

    @Test
    public void testParseMdTagMismatchPositionsThrowsOnInvalidCharacter()
    {
        assertThrows(IllegalArgumentException.class, () -> parseMdTagMismatchPositions(100, "10*5"));
    }
}
