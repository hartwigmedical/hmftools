package com.hartwig.hmftools.esvee.common.saga;

import static java.lang.String.format;

import com.hartwig.hmftools.common.genome.region.Orientation;
import com.hartwig.hmftools.common.region.BasePosition;

import org.jetbrains.annotations.NotNull;

public record SagaBreakend(
        BasePosition position,
        Orientation orientation
)
{
    public static SagaBreakend fromString(final String string)
    {
        String[] parts = string.split(":");
        if(parts.length != 3)
        {
            throw new IllegalArgumentException("Invalid breakend string: " + string);
        }
        String chromosome = parts[0];
        int position = Integer.parseInt(parts[1]);
        Orientation orientation = Orientation.fromByteStr(parts[2]);
        return new SagaBreakend(new BasePosition(chromosome, position), orientation);
    }

    @Override
    @NotNull
    public String toString()
    {
        return format("%s:%s", position, orientation);
    }
}
