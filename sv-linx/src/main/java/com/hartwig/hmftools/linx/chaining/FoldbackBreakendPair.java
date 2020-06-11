package com.hartwig.hmftools.linx.chaining;

import static com.hartwig.hmftools.linx.analysis.SvUtilities.formatJcn;

import java.util.List;

import com.hartwig.hmftools.linx.types.SvBreakend;

public class FoldbackBreakendPair
{
    public final SvBreakend BreakendStart;
    public final SvBreakend BreakendEnd;
    public double Ploidy;
    public final SvChain Chain;

    public FoldbackBreakendPair(final SvBreakend breakend1, final SvBreakend breakend2, double ploidy, final SvChain chain)
    {
        BreakendStart = breakend1;
        BreakendEnd = breakend2;
        Ploidy = ploidy;
        Chain = chain;
    }

    public boolean isChained()
    {
        return Chain != null;
    }

    public String toString()
    {
        return String.format("%s & %s ploidy(%s)", BreakendStart, BreakendEnd, formatJcn(Ploidy));
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

    public static void updateBreakendPair(final List<FoldbackBreakendPair> pairs, final FoldbackBreakendPair pair)
    {
        for(int i = 0; i < pairs.size(); ++i)
        {
            if(pairs.get(i).matches(pair))
            {
                pairs.set(i, pair);
                return;
            }
        }
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
