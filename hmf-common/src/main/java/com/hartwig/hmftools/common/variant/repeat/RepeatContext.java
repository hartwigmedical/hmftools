package com.hartwig.hmftools.common.variant.repeat;

import org.apache.logging.log4j.util.Strings;
import org.jetbrains.annotations.NotNull;

public class RepeatContext
{
    private final byte[] bases;
    private final int repeatIndex;
    private final int startIndex;
    private final int endIndex;
    private final int length;
    private final int forwardCount;
    private final int backwardCount;

    public RepeatContext(final byte[] bases, final int repeatIndex, final int startIndex, final int endIndex, final int length,
            int forwardCount, int backwardCount)
    {
        this.bases = bases;
        this.repeatIndex = repeatIndex;
        this.length = length;
        this.backwardCount = backwardCount;
        this.forwardCount = forwardCount;
        this.startIndex = startIndex;
        this.endIndex = endIndex;
    }

    public int startIndex()
    {
        return startIndex;
    }

    public int endIndex()
    {
        return endIndex;
    }

    public int count()
    {
        return forwardCount + backwardCount;
    }

    @NotNull
    public String sequence()
    {
        return toString();
    }

    @Override
    public String toString()
    {
        return length > 0 ? new String(bases, repeatIndex, length) : Strings.EMPTY;
    }
}