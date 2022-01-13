package com.hartwig.hmftools.sage.read;

import static java.util.Arrays.fill;

import static com.hartwig.hmftools.sage.SageConstants.MATCHING_BASE_QUALITY;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertFalse;
import static org.junit.Assert.assertTrue;

import com.hartwig.hmftools.sage.common.IndexedBases;

import org.apache.logging.log4j.util.Strings;
import org.jetbrains.annotations.NotNull;
import org.junit.Test;

public class ReadContextTest
{
    private ReadContext expandClonedReadContext(final ReadContext readContext, int leftCentreIndex, int rightCentreIndex)
    {
        ReadContext newContext = new ReadContext(
                readContext.Position, readContext.Repeat, readContext.RepeatCount, readContext.Microhomology,
                readContext.indexedBases(), readContext.baseQualities(), readContext.hasIncompleteCore());

        newContext.extendCore(leftCentreIndex, rightCentreIndex);
        return newContext;
    }

    @Test
    public void testExpand()
    {
        final ReadContext initial = create(1000, 2, "CAT", "ACGCA", "GT", 2);

        assertEquals("ACGCA", initial.centerBases());
        assertEquals("TACGCA", expandClonedReadContext(
                initial, initial.readBasesLeftCentreIndex() - 1, initial.readBasesRightCentreIndex()).centerBases());

        assertEquals("ACGCAG", expandClonedReadContext(
                initial, initial.readBasesLeftCentreIndex(), initial.readBasesRightCentreIndex() + 1).centerBases());

        final ReadContext allComplete = expandClonedReadContext(
                initial, initial.readBasesLeftCentreIndex() - 3, initial.readBasesRightCentreIndex() + 2);

        assertFalse(allComplete.hasIncompleteCore());
        assertEquals("CATACGCAGT", allComplete.centerBases());

        final ReadContext allIncompleteLeft = expandClonedReadContext(
                initial, initial.readBasesLeftCentreIndex() - 4, initial.readBasesRightCentreIndex() + 2);

        assertTrue(allIncompleteLeft.hasIncompleteCore());
        assertEquals("CATACGCAGT", allIncompleteLeft.centerBases());

        final ReadContext allIncompleteRight =
                expandClonedReadContext(
                        initial, initial.readBasesLeftCentreIndex() - 3, initial.readBasesRightCentreIndex() + 3);

        assertTrue(allIncompleteRight.hasIncompleteCore());
        assertEquals("CATACGCAGT", allIncompleteRight.centerBases());
    }

    @NotNull
    public static ReadContext simpleSnv(int refPosition, final String leftFlank, final String core,
            final String rightFlank)
    {
        assert (core.length() == 5);
        return create(refPosition, 2, leftFlank, core, rightFlank);
    }

    public static ReadContext create(int refPosition, int indexInCore, final String leftFlank, final String core, final String rightFlank)
    {
        return create(refPosition, indexInCore, leftFlank, core, rightFlank, Math.max(leftFlank.length(), rightFlank.length()));
    }

    public static ReadContext create(
            int refPosition, int indexInCore, final String leftFlank, final String core, final String rightFlank, int flankSize)
    {
        if(indexInCore >= core.length())
        {
            throw new IllegalArgumentException();
        }

        int leftCentreIndex = leftFlank.length();
        int rightCentreIndex = leftCentreIndex + core.length() - 1;
        int readIndex = leftFlank.length() - 1 + indexInCore;
        String readBases = leftFlank + core + rightFlank;
        int coreFlankLength = readBases.length();

        /*
        return new ReadContext("",
                refPosition,
                readIndex,
                leftFlank.length(),
                rightCentreIndex,
                flankSize,
                (leftFlank + core + rightFlank).getBytes(),
                Strings.EMPTY);
        */

        int adjLeftCentreIndex = Math.max(leftCentreIndex, 0);
        int adjRightCentreIndex = Math.min(rightCentreIndex, coreFlankLength - 1);
        boolean incompleteCore = adjLeftCentreIndex != leftCentreIndex || adjRightCentreIndex != rightCentreIndex;

        IndexedBases readBasesIndexed = new IndexedBases(refPosition, readIndex, adjLeftCentreIndex, adjRightCentreIndex, flankSize, readBases.getBytes());

        int[] baseQualities = makeDefaultBaseQualitities(coreFlankLength);

        return new ReadContext(refPosition, "", 0, "", readBasesIndexed, baseQualities, incompleteCore);
    }

    public static int[] makeDefaultBaseQualitities(int coreFlankLength)
    {
        int[] qualities = new int[coreFlankLength];
        fill(qualities, MATCHING_BASE_QUALITY);
        return qualities;

    }
}
