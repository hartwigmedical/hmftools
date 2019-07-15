package com.hartwig.hmftools.linx.chaining;

import java.util.List;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.linx.types.SvLinkedPair;

public class ProposedLinks
{
    public final List<SvLinkedPair> Links;
    public final SvChain ChainTarget;
    public final double Ploidy;

    public final String LinkType;

    public static String PL_TYPE_NONE = "None";
    public static String PL_TYPE_FOLDBACK = "Foldback";
    public static String PL_TYPE_COMPLEX_DUP = "ComplexDup";

    public ProposedLinks(List<SvLinkedPair> links, double ploidy, final SvChain targetChain, final String linkType)
    {
        Links = links;
        Ploidy = ploidy;
        ChainTarget = targetChain;
        LinkType = linkType;
    }

    public ProposedLinks(SvLinkedPair link, double ploidy)
    {
        Links = Lists.newArrayList(link);
        Ploidy = ploidy;
        ChainTarget = null;
        LinkType = PL_TYPE_NONE;
    }

}
