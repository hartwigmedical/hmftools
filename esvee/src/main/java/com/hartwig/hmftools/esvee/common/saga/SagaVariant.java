package com.hartwig.hmftools.esvee.common.saga;

import static java.lang.String.format;

import java.util.stream.Stream;

import org.jetbrains.annotations.NotNull;

public record SagaVariant(
        String id,
        SagaBreakend breakend1,
        SagaBreakend breakend2
)
{
    public Stream<SagaBreakend> breakends()
    {
        return Stream.of(breakend1, breakend2);
    }

    @Override
    @NotNull
    public String toString()
    {
        return format("%s %s-%s", id, breakend1, breakend2);
    }

    public String noSpaceString()
    {
        String result = format("%s:%s-%s", id, breakend1, breakend2);
        result = result.replace(" ", "_");
        return result;
    }
}
