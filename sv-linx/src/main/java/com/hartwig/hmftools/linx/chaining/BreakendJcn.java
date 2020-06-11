package com.hartwig.hmftools.linx.chaining;

import static java.lang.Math.max;

import java.util.List;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.utils.Doubles;

public class BreakendJcn
{
    public final double UnchainedJcn;
    public final double TotalChainJcn;
    public final SvChain MaxJcnChain;
    public final List<SvChain> Chains;
    public final boolean OtherBreakendExhausted;

    public BreakendJcn(double unchainedJcn, boolean otherBreakendExhausted)
    {
        UnchainedJcn = unchainedJcn;
        TotalChainJcn = 0;
        MaxJcnChain = null;
        Chains = Lists.newArrayList();
        OtherBreakendExhausted = otherBreakendExhausted;
    }

    public BreakendJcn(double unchainedJcn, double totalChainJcn, final SvChain maxJcnChain,
            final List<SvChain> chains, boolean otherBreakendExhausted)
    {
        UnchainedJcn = unchainedJcn;
        TotalChainJcn = totalChainJcn;
        MaxJcnChain = maxJcnChain;
        Chains = chains;
        OtherBreakendExhausted = otherBreakendExhausted;
    }

    public double unlinkedJcn()
    {
        if(TotalChainJcn > 0)
        {
            return max(MaxJcnChain.jcn(), UnchainedJcn);
        }
        else
        {
            return UnchainedJcn;
        }
    }

    public boolean multiConnections()
    {
        return Chains.size() > 1 || (!Chains.isEmpty() && UnchainedJcn > 0);
    }

    public boolean exhaustedInChain()
    {
        return Chains.size() == 1 && Doubles.equal(UnchainedJcn, 0);
    }
}
