package com.hartwig.hmftools.linx.chaining;

import static java.lang.Math.min;

import static com.hartwig.hmftools.linx.analysis.SvUtilities.copyNumbersEqual;
import static com.hartwig.hmftools.linx.chaining.ChainPloidyLimits.ploidyOverlap;
import static com.hartwig.hmftools.linx.chaining.ChainingRule.calcRulePriority;

import java.util.List;
import java.util.Map;

import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import com.hartwig.hmftools.linx.types.SvBreakend;
import com.hartwig.hmftools.linx.types.SvLinkedPair;

public class ProposedLinks
{
    public final List<SvLinkedPair> Links;

    private final List<ChainingRule> mRules;
    private int mPriority;

    private SvChain mChainTarget;
    private String mChainConnectType;

    private double mPloidy; // the min ploidy if the breakends differ, otherwise the median
    private Map<SvBreakend, Double> mBreakendMinPloidy;

    private String mPloidyMatchType;
    private long mShortestDistance;

    public static String CONN_TYPE_NONE = "None";
    public static String CONN_TYPE_FOLDBACK = "Foldback";
    public static String CONN_TYPE_COMPLEX_DUP = "ComplexDup";

    public static String PM_MATCHED = "Matched";
    public static String PM_NONE = "None";
    public static String PM_OVERLAP = "Overlap";

    public ProposedLinks(SvLinkedPair link, ChainingRule rule)
    {
        Links = Lists.newArrayList(link);

        mRules = Lists.newArrayList(rule);
        mPriority = calcRulePriority(mRules);

        mPloidy = 0;
        mBreakendMinPloidy = Maps.newHashMap();
        mChainTarget = null;
        mChainConnectType = CONN_TYPE_NONE;
        mPloidyMatchType = PM_NONE;

        mShortestDistance = link.length();
    }

    public ProposedLinks(List<SvLinkedPair> links, ChainingRule rule, final SvChain targetChain, final String connectType)
    {
        Links = links;

        mRules = Lists.newArrayList(rule);
        mPriority = calcRulePriority(mRules);

        mPloidy = 0;
        mBreakendMinPloidy = Maps.newHashMap();
        mChainTarget = targetChain;
        mChainConnectType = connectType;

        mShortestDistance = Links.stream().mapToLong(x -> x.length()).min().getAsLong();
    }

    public final SvChain targetChain() { return mChainTarget; }
    public final String chainConnectType() { return mChainConnectType; };

    public double ploidy() { return mPloidy; }

    public void addBreakendPloidies(final SvBreakend breakend1, double ploidy1, final SvBreakend breakend2, double ploidy2, boolean recalc)
    {
        addBreakendPloidies(breakend1, ploidy1, breakend2, ploidy2);

        if(!recalc)
            return;

        if(copyNumbersEqual(ploidy1, ploidy2))
        {
            mPloidyMatchType = PM_MATCHED;
        }
        else if(ploidyOverlap(ploidy1, breakend1.getSV().ploidyUncertainty(),
                ploidy2, breakend2.getSV().ploidyUncertainty()))
        {
            mPloidyMatchType = PM_OVERLAP;
        }

        if(mPloidyMatchType == PM_NONE)
            mPloidy = min(ploidy1, ploidy2);
        else
            mPloidy = (ploidy1 + ploidy2) * 0.5;
    }

    public void addBreakendPloidies(final SvBreakend breakend1, double ploidy1, final SvBreakend breakend2, double ploidy2)
    {
        mBreakendMinPloidy.put(breakend1, ploidy1);
        mBreakendMinPloidy.put(breakend2, ploidy2);
    }

    public double breakendPloidy(final SvBreakend breakend)
    {
        return mBreakendMinPloidy.get(breakend);
    }

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
