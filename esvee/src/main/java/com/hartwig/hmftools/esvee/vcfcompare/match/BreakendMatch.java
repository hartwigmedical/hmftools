package com.hartwig.hmftools.esvee.vcfcompare.match;

import java.util.StringJoiner;

import com.hartwig.hmftools.esvee.vcfcompare.VariantBreakend;

import org.jetbrains.annotations.Nullable;

public class BreakendMatch
{
    public @Nullable VariantBreakend OldBreakend;
    public @Nullable VariantBreakend NewBreakend;
    public MatchType Type;

    public BreakendMatch(@Nullable VariantBreakend oldBreakend, @Nullable VariantBreakend newBreakend, MatchType type)
    {
        if(oldBreakend == null && newBreakend == null)
            throw new IllegalArgumentException("`oldBreakend` and `newBreakend` cannot both be null");

        OldBreakend = oldBreakend;
        NewBreakend = newBreakend;
        Type = type;
    }

    @Override
    public String toString()
    {
        StringJoiner sj = new StringJoiner(" ");

        if(OldBreakend != null)
            sj.add(String.format("oldBreakend(%s)", OldBreakend.coordStr()));

        if(NewBreakend != null)
            sj.add(String.format("newBreakend(%s)", NewBreakend.coordStr()));

        return sj.toString();
    }
}
