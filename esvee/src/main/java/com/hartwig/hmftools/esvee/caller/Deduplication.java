package com.hartwig.hmftools.esvee.caller;

import java.util.List;
import java.util.Map;

import com.hartwig.hmftools.esvee.common.FilterType;

public final class Deduplication
{
    public static void deduplicateVariants(final Map<String,List<Breakend>> chrBreakendMap)
    {
        chrBreakendMap.values().forEach(x -> deduplicateBreakends(x));
    }

    private static void deduplicateBreakends(final List<Breakend> breakends)
    {
        for(int i = 0; i < breakends.size() - 1; ++i)
        {
            Breakend breakend = breakends.get(i);

            if(breakend.sv().isFiltered())
                continue;

            for(int j = i + 1; j < breakends.size(); ++j)
            {
                Breakend nextBreakend = breakends.get(j);

                if(nextBreakend.sv().isFiltered()) // skip setting the duplicate filter either are filtered
                    continue;

                if(!breakend.Chromosome.equals(nextBreakend.Chromosome) || nextBreakend.Position > breakend.Position)
                    break;

                if(isDuplicate(breakend, nextBreakend))
                {
                    // keep passing breakend with most support
                    if(keepFirst(breakend, nextBreakend))
                        nextBreakend.sv().addFilter(FilterType.DUPLICATE);
                    else
                        breakend.sv().addFilter(FilterType.DUPLICATE);
                }
            }
        }
    }

    private static boolean isDuplicate(final Breakend first, final Breakend second)
    {
        if(first.Position != second.Position || first.Orient != second.Orient || !first.InsertSequence.equals(second.InsertSequence))
            return false;

        if(first.isSgl() != second.isSgl())
            return false;

        if(first.isSgl())
            return true;

        Breakend firstOther = first.otherBreakend();
        Breakend secondOther = second.otherBreakend();

        return firstOther.Position == secondOther.Position
            && firstOther.Orient == secondOther.Orient
            && firstOther.Chromosome.equals(secondOther.Chromosome);
    }

    private static boolean keepFirst(final Breakend first, final Breakend second)
    {
        int firstSupport = first.fragmentCount();
        int secondSupport = second.fragmentCount();

        if(firstSupport != secondSupport)
            return firstSupport > secondSupport;

        return first.sv().qual() >= second.sv().qual();
    }
}
