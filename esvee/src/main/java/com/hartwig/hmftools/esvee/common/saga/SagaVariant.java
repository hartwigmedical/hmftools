package com.hartwig.hmftools.esvee.common.saga;

import static java.lang.String.format;

import java.util.stream.Stream;

import org.jetbrains.annotations.NotNull;

public record SagaVariant(
        String id,
        SagaBreakend breakend1,
        SagaBreakend breakend2,
        String insertSequence
)
{
    public SagaVariant
    {
        // Only indels should be in the SAGA resource, and downstream code assumes all SAGA variants are indels.
        if(!breakend1.chromosome().equals(breakend2.chromosome()))
        {
            throw new IllegalArgumentException();
        }
        if(!(breakend1.orientation().isForward() && breakend2.orientation().isReverse()))
        {
            throw new IllegalArgumentException();
        }
        if(!(breakend1.position() < breakend2.position()))
        {
            throw new IllegalArgumentException();
        }
        boolean isInsert = breakend1.position() + 1 == breakend2.position();
        if(isInsert && insertSequence.isEmpty())
        {
            throw new IllegalArgumentException();
        }
    }

    public Stream<SagaBreakend> breakends()
    {
        return Stream.of(breakend1, breakend2);
    }

    public boolean isInsert()
    {
        return !insertSequence().isEmpty();
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
