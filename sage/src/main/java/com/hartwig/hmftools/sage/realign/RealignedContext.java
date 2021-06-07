package com.hartwig.hmftools.sage.realign;

import org.jetbrains.annotations.NotNull;

public class RealignedContext
{

    private final RealignedType type;
    private final int repeatCount;

    public RealignedContext(@NotNull final RealignedType type, final int repeatCount)
    {
        this.type = type;
        this.repeatCount = repeatCount;
    }

    @NotNull
    public RealignedType type()
    {
        return type;
    }

    public int repeatCount()
    {
        return repeatCount;
    }
}
