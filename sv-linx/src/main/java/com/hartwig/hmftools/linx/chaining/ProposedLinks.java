package com.hartwig.hmftools.linx.chaining;

import java.util.List;

import com.hartwig.hmftools.linx.types.SvLinkedPair;

public class ProposedLinks
{
    public final List<SvLinkedPair> Links;
    public final SvChain ChainTarget;
    public final double Ploidy;

    public ProposedLinks(List<SvLinkedPair> links, double ploidy, final SvChain targetChain)
    {
        Links = links;
        Ploidy = ploidy;
        ChainTarget = targetChain;
    }

}
