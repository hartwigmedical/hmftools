package com.hartwig.hmftools.sage.common;

import static java.util.Arrays.fill;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertFalse;
import static org.junit.Assert.assertTrue;

import com.hartwig.hmftools.sage.candidate_.IndexedBases_;
import com.hartwig.hmftools.sage.candidate_.ReadContext_;

import org.jetbrains.annotations.NotNull;
import org.junit.Test;

public class ReadContextTest
{
    private ReadContext_ expandClonedReadContext(final ReadContext_ readContext, int leftCentreIndex, int rightCentreIndex)
    {
        ReadContext_ newContext = new ReadContext_(
                readContext.Position, readContext.Repeat, readContext.RepeatCount, readContext.Microhomology,
                readContext.indexedBases(), readContext.hasIncompleteCore());

        newContext.extendCore(leftCentreIndex, rightCentreIndex);
        return newContext;
    }

    @Test
    public void testExpand()
    {
        final ReadContext_ initial = create(1000, 2, "CAT", "ACGCA", "GT", 2);

        assertEquals("ACGCA", initial.coreString());
        assertEquals("TACGCA", expandClonedReadContext(
                initial, initial.readBasesLeftCentreIndex() - 1, initial.readBasesRightCentreIndex()).coreString());

        assertEquals("ACGCAG", expandClonedReadContext(
                initial, initial.readBasesLeftCentreIndex(), initial.readBasesRightCentreIndex() + 1).coreString());

        final ReadContext_ allComplete = expandClonedReadContext(
                initial, initial.readBasesLeftCentreIndex() - 3, initial.readBasesRightCentreIndex() + 2);

        assertFalse(allComplete.hasIncompleteCore());
        assertEquals("CATACGCAGT", allComplete.coreString());

        final ReadContext_ allIncompleteLeft = expandClonedReadContext(
                initial, initial.readBasesLeftCentreIndex() - 4, initial.readBasesRightCentreIndex() + 2);

        assertTrue(allIncompleteLeft.hasIncompleteCore());
        assertEquals("CATACGCAGT", allIncompleteLeft.coreString());

        final ReadContext_ allIncompleteRight =
                expandClonedReadContext(
                        initial, initial.readBasesLeftCentreIndex() - 3, initial.readBasesRightCentreIndex() + 3);

        assertTrue(allIncompleteRight.hasIncompleteCore());
        assertEquals("CATACGCAGT", allIncompleteRight.coreString());
    }

    @NotNull
    public static ReadContext_ simpleSnv(int refPosition, final String leftFlank, final String core,
            final String rightFlank)
    {
        assert (core.length() == 5);
        return create(refPosition, 2, leftFlank, core, rightFlank);
    }

    public static ReadContext_ create(int refPosition, int indexInCore, final String leftFlank, final String core, final String rightFlank)
    {
        return create(refPosition, indexInCore, leftFlank, core, rightFlank, Math.max(leftFlank.length(), rightFlank.length()));
    }

    public static ReadContext_ create(
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

        int adjLeftCentreIndex = Math.max(leftCentreIndex, 0);
        int adjRightCentreIndex = Math.min(rightCentreIndex, coreFlankLength - 1);
        boolean incompleteCore = adjLeftCentreIndex != leftCentreIndex || adjRightCentreIndex != rightCentreIndex;

        IndexedBases_ readBasesIndexed = new IndexedBases_(refPosition, readIndex, adjLeftCentreIndex, adjRightCentreIndex, flankSize, readBases.getBytes());

        return new ReadContext_(refPosition, "", 0, "", readBasesIndexed, incompleteCore);
    }
}
