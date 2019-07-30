package com.hartwig.hmftools.linx.chaining;

import java.util.List;

import com.hartwig.hmftools.linx.types.SvBreakend;

public class FoldbackBreakendPair
{
    public final SvBreakend BreakendStart;
    public final SvBreakend BreakendEnd;
    public final double Ploidy;

    public FoldbackBreakendPair(final SvBreakend breakend1, final SvBreakend breakend2, double ploidy)
    {
        BreakendStart = breakend1;
        BreakendEnd = breakend2;
        Ploidy = ploidy;
    }

    public boolean isChained()
    {
        return BreakendStart.getSV() != BreakendEnd.getSV() || BreakendStart == BreakendEnd;
    }

    public String toString()
    {
        return String.format("%s & %s ploidy(%.1f)", BreakendStart, BreakendEnd, Ploidy);
    }

    public boolean matches(final FoldbackBreakendPair other)
    {
        return (other.BreakendEnd == BreakendEnd && other.BreakendStart == BreakendStart)
                || (other.BreakendStart == BreakendEnd && other.BreakendEnd == BreakendStart);
    }

    public static void addByPloidy(final List<FoldbackBreakendPair> pairs, FoldbackBreakendPair newPair)
    {
        int index = 0;
        while(index < pairs.size())
        {
            if(pairs.get(index).Ploidy < newPair.Ploidy)
                break;

            ++index;
        }

        pairs.add(index, newPair);
    }

    public static boolean containsBreakendPair(final List<FoldbackBreakendPair> pairs, final FoldbackBreakendPair pair)
    {
        return pairs.stream().anyMatch(x -> x.matches(pair));
    }

    public static void removeBreakendPair(final List<FoldbackBreakendPair> pairs, final FoldbackBreakendPair pair)
    {
        for(int i = 0; i < pairs.size(); ++i)
        {
            if(pairs.get(i).matches(pair))
            {
                pairs.remove(i);
                return;
            }
        }
    }
}
