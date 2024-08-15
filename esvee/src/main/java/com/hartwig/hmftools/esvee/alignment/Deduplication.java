package com.hartwig.hmftools.esvee.alignment;

import java.util.List;

import com.hartwig.hmftools.esvee.common.FilterType;

public final class Deduplication
{
    public static void deduplicateBreakends(final List<Breakend> breakends)
    {
        for(int i = 0; i < breakends.size() - 1; ++i)
        {
            Breakend breakend = breakends.get(i);

            if(!breakend.passing())
                continue;

            for(int j = i + 1; j < breakends.size(); ++j)
            {
                Breakend nextBreakend = breakends.get(j);

                if(!breakend.Chromosome.equals(nextBreakend.Chromosome) || nextBreakend.Position > breakend.Position)
                    break;

                if(isDuplicate(breakend, nextBreakend))
                {
                    // keep breakend with most support
                    if(keepFirst(breakend, nextBreakend))
                        nextBreakend.addFilter(FilterType.DUPLICATE);
                    else
                        breakend.addFilter(FilterType.DUPLICATE);
                }
            }
        }
    }

    private static boolean isDuplicate(final Breakend first, final Breakend second)
    {
        return first.Orient == second.Orient && first.InsertedBases.equals(second.InsertedBases);
    }

    private static boolean keepFirst(final Breakend first, final Breakend second)
    {
        if(first.isSingle() != second.isSingle())
            return !first.isSingle();

        int firstSupport = first.sampleSupport().stream().mapToInt(x -> x.totalSupport()).sum();
        int secondSupport = second.sampleSupport().stream().mapToInt(x -> x.totalSupport()).sum();

        if(firstSupport != secondSupport)
            return firstSupport > secondSupport;

        return first.calcSvQual() >= second.calcSvQual();
    }
}
