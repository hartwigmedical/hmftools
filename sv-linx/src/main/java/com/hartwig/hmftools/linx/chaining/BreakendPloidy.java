package com.hartwig.hmftools.linx.chaining;

import static java.lang.Math.max;

import java.util.List;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.numeric.Doubles;

public class BreakendPloidy
{
    public final double UnchainedPloidy;
    public final double TotalChainPloidy;
    public final SvChain MaxPloidyChain;
    public final List<SvChain> Chains;
    public final boolean OtherBreakendExhausted;

    public BreakendPloidy(double unchainedPloidy, boolean otherBreakendExhausted)
    {
        UnchainedPloidy = unchainedPloidy;
        TotalChainPloidy = 0;
        MaxPloidyChain = null;
        Chains = Lists.newArrayList();
        OtherBreakendExhausted = otherBreakendExhausted;
    }

    public BreakendPloidy(double unchainedPloidy, double totalChainPloidy, final SvChain maxPloidyChain,
            final List<SvChain> chains, boolean otherBreakendExhausted)
    {
        UnchainedPloidy = unchainedPloidy;
        TotalChainPloidy = totalChainPloidy;
        MaxPloidyChain = maxPloidyChain;
        Chains = chains;
        OtherBreakendExhausted = otherBreakendExhausted;
    }

    public double unlinkedPloidy()
    {
        if(TotalChainPloidy > 0)
        {
            return max(MaxPloidyChain.ploidy(), UnchainedPloidy);
        }
        else
        {
            return UnchainedPloidy;
        }
    }

    public boolean multiConnections()
    {
        return Chains.size() > 1 || (!Chains.isEmpty() && UnchainedPloidy > 0);
    }

    public boolean exhaustedInChain()
    {
        return Chains.size() == 1 && Doubles.equal(UnchainedPloidy, 0);
    }
}
