package com.hartwig.hmftools.linx.chaining;

import static java.lang.Math.abs;

import static com.hartwig.hmftools.linx.chaining.ChainingRule.calcRulePriority;

import java.util.List;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.linx.types.SvLinkedPair;

public class ProposedLinks
{
    public final List<SvLinkedPair> Links;
    public final double Ploidy;

    private final List<ChainingRule> mRules;
    private int mPriority;

    private SvChain mChainTarget;
    private String mChainConnectType;
    private String mPloidyMatchType;
    private long mShortestDistance;

    public static String PL_TYPE_NONE = "None";
    public static String PL_TYPE_FOLDBACK = "Foldback";
    public static String PL_TYPE_COMPLEX_DUP = "ComplexDup";

    public static String PM_MATCHED = "Matched";
    public static String PM_NONE = "None";
    public static String PM_OVERLAP = "Overlap";

    public ProposedLinks(SvLinkedPair link, double ploidy, ChainingRule rule)
    {
        Links = Lists.newArrayList(link);
        Ploidy = ploidy;
        mRules = Lists.newArrayList(rule);
        mPriority = calcRulePriority(mRules);
        mChainTarget = null;
        mChainConnectType = PL_TYPE_NONE;
        mPloidyMatchType = PM_NONE;

        mShortestDistance = link.length();
    }

    public ProposedLinks(List<SvLinkedPair> links, double ploidy, ChainingRule rule, final SvChain targetChain, final String connectType)
    {
        Links = links;
        Ploidy = ploidy;
        mRules = Lists.newArrayList(rule);
        mPriority = calcRulePriority(mRules);

        mChainTarget = targetChain;
        mChainConnectType = connectType;

        mShortestDistance = Links.stream().mapToLong(x -> x.length()).min().getAsLong();
    }

    public final SvChain targetChain() { return mChainTarget; }
    public final String chainConnectType() { return mChainConnectType; };

    public void setPloidyMatch(final String type) { mPloidyMatchType = type; }
    public final String ploidyMatchType() { return mPloidyMatchType; }

    public void addRule(ChainingRule rule)
    {
        if(mRules.contains(rule))
            return;

        mRules.add(rule);
        mPriority = calcRulePriority(mRules);
    }

    public int priority() { return mPriority; }

    public ChainingRule topRule()
    {
        return mRules.get(0);
    }

    public long shortestLinkDistance() { return mShortestDistance; }

}
