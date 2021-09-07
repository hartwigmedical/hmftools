package com.hartwig.hmftools.sage.read;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertFalse;
import static org.junit.Assert.assertTrue;

import org.apache.logging.log4j.util.Strings;
import org.jetbrains.annotations.NotNull;
import org.junit.Test;

public class ReadContextTest
{

    @Test
    public void testExpand()
    {
        final ReadContext initial = create(1000, 2, "CAT", "ACGCA", "GT", 2);
        assertEquals("ACGCA", initial.centerBases());
        assertEquals("TACGCA", initial.extend(initial.readBasesLeftCentreIndex() - 1, initial.readBasesRightCentreIndex()).centerBases());
        assertEquals("ACGCAG", initial.extend(initial.readBasesLeftCentreIndex(), initial.readBasesRightCentreIndex() + 1).centerBases());

        final ReadContext allComplete = initial.extend(initial.readBasesLeftCentreIndex() - 3, initial.readBasesRightCentreIndex() + 2);
        assertFalse(allComplete.incompleteCore());
        assertEquals("CATACGCAGT", allComplete.centerBases());

        final ReadContext allIncompleteLeft =
                initial.extend(initial.readBasesLeftCentreIndex() - 4, initial.readBasesRightCentreIndex() + 2);
        assertTrue(allIncompleteLeft.incompleteCore());
        assertEquals("CATACGCAGT", allIncompleteLeft.centerBases());

        final ReadContext allIncompleteRight =
                initial.extend(initial.readBasesLeftCentreIndex() - 3, initial.readBasesRightCentreIndex() + 3);
        assertTrue(allIncompleteRight.incompleteCore());
        assertEquals("CATACGCAGT", allIncompleteRight.centerBases());
    }

    @NotNull
    public static ReadContext simpleSnv(int refPosition, @NotNull final String leftFlank, @NotNull final String core,
            @NotNull final String rightFlank)
    {
        assert (core.length() == 5);
        return create(refPosition, 2, leftFlank, core, rightFlank);
    }

    @NotNull
    public static ReadContext create(int refPosition, int indexInCore, @NotNull final String leftFlank, @NotNull final String core,
            @NotNull final String rightFlank)
    {
        return create(refPosition, indexInCore, leftFlank, core, rightFlank, Math.max(leftFlank.length(), rightFlank.length()));
    }

    @NotNull
    public static ReadContext create(int refPosition, int indexInCore, @NotNull final String leftFlank, @NotNull final String core,
            @NotNull final String rightFlank, int flankSize)
    {
        if(indexInCore >= core.length())
        {
            throw new IllegalArgumentException();
        }

        int leftCentreIndex = leftFlank.length();
        int rightCentreIndex = leftCentreIndex + core.length() - 1;
        int index = leftFlank.length() - 1 + indexInCore;

        return new ReadContext("",
                refPosition,
                index,
                leftFlank.length(),
                rightCentreIndex,
                flankSize,
                (leftFlank + core + rightFlank).getBytes(),
                Strings.EMPTY);
    }
}
