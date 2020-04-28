package com.hartwig.hmftools.sage.read;

import org.apache.logging.log4j.util.Strings;
import org.jetbrains.annotations.NotNull;

public class ReadContextTest {

    @NotNull
    public static ReadContext simpleSnv(int refPosition, @NotNull final String leftFlank, @NotNull final String core,
            @NotNull final String rightFlank) {
        assert (core.length() == 5);
        return create(refPosition, 2, leftFlank, core, rightFlank);
    }

    @NotNull
    public static ReadContext create(int refPosition, int indexInCore, @NotNull final String leftFlank, @NotNull final String core,
            @NotNull final String rightFlank) {
        if (indexInCore >= core.length()) {
            throw new IllegalArgumentException();
        }

        int leftCentreIndex = leftFlank.length();
        int rightCentreIndex = leftCentreIndex + core.length() - 1;
        int flankSize = Math.max(leftFlank.length(), rightFlank.length());
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
