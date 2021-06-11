package com.hartwig.hmftools.sage.realign;

import org.jetbrains.annotations.NotNull;

public class RealignedContext
{
    public final RealignedType Type;
    public final int RepeatCount;

    public RealignedContext(@NotNull final RealignedType type, final int repeatCount)
    {
        Type = type;
        RepeatCount = repeatCount;
    }
}
