package com.hartwig.hmftools.linx.chaining;

import static com.hartwig.hmftools.linx.analysis.SvUtilities.copyNumbersEqual;
import static com.hartwig.hmftools.linx.analysis.SvUtilities.formatPloidy;
import static com.hartwig.hmftools.linx.chaining.ChainPloidyLimits.ploidyOverlap;
import static com.hartwig.hmftools.linx.chaining.ChainingRule.COMP_DUP_SPLIT;
import static com.hartwig.hmftools.linx.chaining.ChainingRule.FOLDBACK_SPLIT;
import static com.hartwig.hmftools.linx.chaining.ChainingRule.calcRulePriority;

import java.util.List;
import java.util.Map;

import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import com.hartwig.hmftools.common.utils.Doubles;
import com.hartwig.hmftools.linx.types.SvBreakend;
import com.hartwig.hmftools.linx.types.SvLinkedPair;

public class ProposedLinks
{
    public final List<SvLinkedPair> Links;

    private final List<ChainingRule> mRules;
    private int mPriority;

    // for complex types, the chain which is being replicated by the foldback or complex dup
    private SvChain mChainTarget;

    // chain operating as a foldback - only relevant for foldback-splitting links
    private SvChain mFoldbackChain;

    private double mPloidy; // the min ploidy if the breakends differ, otherwise the average
    private Map<SvBreakend, Double> mBreakendPloidy;

    // indication of whether a breakend is exhausted by making these links
    private Map<SvBreakend, Boolean> mExhaustBreakend;

    private String mPloidyMatchType;
    private int mShortestDistance;

    public static String PM_MATCHED = "Matched";
    public static String PM_NONE = "None";
    public static String PM_OVERLAP = "Overlap";

    public ProposedLinks(SvLinkedPair link, ChainingRule rule)
    {
        Links = Lists.newArrayList(link);

        mRules = Lists.newArrayList(rule);
        mPriority = calcRulePriority(mRules);

        mPloidy = 0;
        mBreakendPloidy = Maps.newHashMap();
        mExhaustBreakend = Maps.newHashMap();
        mChainTarget = null;
        mFoldbackChain = null;
        mPloidyMatchType = PM_NONE;

        mShortestDistance = link.length();
    }

    public ProposedLinks(List<SvLinkedPair> links, ChainingRule rule, final SvChain targetChain, final SvChain foldbackChain)
    {
        Links = links;

        mRules = Lists.newArrayList(rule);
        mPriority = calcRulePriority(mRules);

        mPloidy = 0;
        mBreakendPloidy = Maps.newHashMap();
        mExhaustBreakend = Maps.newHashMap();
        mChainTarget = targetChain;
        mFoldbackChain = foldbackChain;
        mPloidyMatchType = PM_NONE;

        mShortestDistance = Links.stream().mapToInt(x -> x.length()).min().getAsInt();
    }

    public final SvChain targetChain() { return mChainTarget; }
    public final SvChain foldbackChain() { return mFoldbackChain; }

    public final ChainingRule getSplittingRule()
    {
        if(mRules.contains(FOLDBACK_SPLIT))
            return FOLDBACK_SPLIT;
        else
            return COMP_DUP_SPLIT;
    }

    public final boolean multiConnection() { return Links.size() > 1; };

    public double ploidy() { return mPloidy; }

    public void setLowerPloidy(double ploidy)
    {
        if(ploidy >= mPloidy || multiConnection())
            return;

        mPloidy = ploidy;

        // this may change whether a link is exhausted
        for(Map.Entry<SvBreakend,Double> entry : mBreakendPloidy.entrySet())
        {
            final SvBreakend breakend = entry.getKey();
            double bePloidy = entry.getValue();
            mExhaustBreakend.put(breakend, Doubles.equal(bePloidy, mPloidy));
        }
    }

    public void addComDupBreakends(
            final SvBreakend compDupStart, final SvBreakend compDupEnd, double compDupPloidy,
            final SvBreakend otherStart, final SvBreakend otherEnd, double otherPloidy)
    {
        if(compDupPloidy == 0 || otherPloidy == 0)
            return;

        // no match for this type of connection

        mPloidy = compDupPloidy;

        mBreakendPloidy.put(compDupStart, compDupPloidy);
        mBreakendPloidy.put(compDupEnd, compDupPloidy);
        mExhaustBreakend.put(compDupStart, true);
        mExhaustBreakend.put(compDupEnd, true);

        mBreakendPloidy.put(otherStart, otherPloidy);
        mBreakendPloidy.put(otherEnd, otherPloidy);
    }

    public void addFoldbackBreakends(
            final SvBreakend foldbackStart, final SvBreakend foldbackEnd, double foldbackPloidy,
            final SvBreakend otherBreakend, double otherPloidy, double otherUncertainty)
    {
        if(foldbackPloidy == 0 || otherPloidy == 0)
            return;

        if (copyNumbersEqual(foldbackPloidy * 2, otherPloidy))
        {
            mPloidyMatchType = PM_MATCHED;
        }
        else if (ploidyOverlap(foldbackPloidy * 2, foldbackStart.ploidyUncertainty(), otherPloidy, otherUncertainty))
        {
            mPloidyMatchType = PM_OVERLAP;
        }

        if(mPloidyMatchType == PM_NONE)
            mPloidy = foldbackPloidy;
        else
            mPloidy = (foldbackPloidy + otherPloidy/2) * 0.5;

        mBreakendPloidy.put(foldbackStart, foldbackPloidy);
        mBreakendPloidy.put(foldbackEnd, foldbackPloidy);

        if(foldbackStart.getSV() == foldbackEnd.getSV() && copyNumbersEqual(foldbackStart.ploidy(), foldbackPloidy))
        {
            // mark the foldback breakends as exhausted only if clear that they are fully used in this pair of links
            mExhaustBreakend.put(foldbackStart, true);
            mExhaustBreakend.put(foldbackEnd, true);
        }

        mBreakendPloidy.put(otherBreakend, otherPloidy);
    }

    public void addBreakendPloidies(
            final SvBreakend breakend1, double ploidy1,
            final SvBreakend breakend2, double ploidy2)
    {
        if(ploidy1 == 0 || ploidy2 == 0)
            return;

        mBreakendPloidy.put(breakend1, ploidy1);
        mBreakendPloidy.put(breakend2, ploidy2);

        if (copyNumbersEqual(ploidy1, ploidy2))
        {
            mPloidyMatchType = PM_MATCHED;
        }
        else if (ploidyOverlap(ploidy1, breakend1.ploidyUncertainty(), ploidy2, breakend2.ploidyUncertainty()))
        {
            mPloidyMatchType = PM_OVERLAP;
        }

        if(mPloidyMatchType == PM_NONE)
        {
            if(ploidy1 > ploidy2)
            {
                mPloidy = ploidy2;
            }
            else
            {
                mPloidy = ploidy1;
            }
        }
        else
        {
            mPloidy = (ploidy1 + ploidy2) * 0.5;
        }
    }

    public double breakendPloidy(final SvBreakend breakend)
    {
        return mBreakendPloidy.get(breakend);
    }

    public boolean exhaustBreakend(final SvBreakend breakend)
    {
        Boolean exhausted = mExhaustBreakend.get(breakend);
        return exhausted != null && exhausted;
    }

    public void overrideBreakendPloidyMatched(final SvBreakend breakend, boolean matched)
    {
        mExhaustBreakend.put(breakend, matched);
    }

    public final String ploidyMatchType() { return mPloidyMatchType; }
    public boolean linkPloidyMatch() { return mPloidyMatchType != PM_NONE; }

    public void addRule(ChainingRule rule)
    {
        if(mRules.contains(rule))
            return;

        mRules.add(rule);
        mPriority = calcRulePriority(mRules);
    }

    public int priority() { return mPriority; }

    public ChainingRule topRule() { return mRules.get(0); }
    public boolean hasRule(ChainingRule rule) { return mRules.contains(rule); }

    public int shortestLinkDistance() { return mShortestDistance; }

    public String toString()
    {
        return String.format("pair(%s) rule(%s) ploidy(%s) match(%s) length(%d) priority(%d)",
                Links.get(0).toString(), mRules.get(0), formatPloidy(mPloidy),
                mPloidyMatchType, mShortestDistance, mPriority);
    }

    public boolean isValid()
    {
        if(Links.isEmpty() || Links.size() > 2)
            return false;

        if(mBreakendPloidy.isEmpty())
            return false;

        if(Links.size() == 1 && mBreakendPloidy.size() != 2)
            return false;

        if(Links.size() == 2 && mBreakendPloidy.size() < 3)
        {
            if(!Links.get(0).matches(Links.get(1)))
                return false;
        }

        return mPloidy > 0;
    }

}
