package com.hartwig.hmftools.linx.chaining;

import java.util.List;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.linx.types.SvLinkedPair;

public class ProposedLinks
{
    public final List<SvLinkedPair> Links;
    public final SvChain ChainTarget;
    public final double Ploidy;
    public final boolean ComplexType;

    public ProposedLinks(List<SvLinkedPair> links, double ploidy, final SvChain targetChain, boolean complexType)
    {
        Links = links;
        Ploidy = ploidy;
        ChainTarget = targetChain;
        ComplexType = complexType;
    }

    public ProposedLinks(SvLinkedPair link, double ploidy)
    {
        Links = Lists.newArrayList(link);
        Ploidy = ploidy;
        ChainTarget = null;
        ComplexType = false;
    }

}
