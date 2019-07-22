package com.hartwig.hmftools.linx.chaining;

import static java.lang.Math.min;

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
import com.hartwig.hmftools.common.numeric.Doubles;
import com.hartwig.hmftools.linx.types.SvBreakend;
import com.hartwig.hmftools.linx.types.SvLinkedPair;

public class ProposedLinks
{
    public final List<SvLinkedPair> Links;

    private final List<ChainingRule> mRules;
    private int mPriority;

    // for complex types, the chain which is being replicated by the foldback or complex dup
    private SvChain mChainTarget;

    private double mPloidy; // the min ploidy if the breakends differ, otherwise the median
    private Map<SvBreakend, Double> mBreakendPloidy;
    private Map<SvBreakend, Boolean> mBreakendPloidyMatched;

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
        mBreakendPloidy = Maps.newHashMap();
        mBreakendPloidyMatched = Maps.newHashMap();
        mChainTarget = null;
        mPloidyMatchType = PM_NONE;

        mShortestDistance = link.length();
    }

    public ProposedLinks(List<SvLinkedPair> links, ChainingRule rule, final SvChain targetChain)
    {
        Links = links;

        mRules = Lists.newArrayList(rule);
        mPriority = calcRulePriority(mRules);

        mPloidy = 0;
        mBreakendPloidy = Maps.newHashMap();
        mBreakendPloidyMatched = Maps.newHashMap();
        mChainTarget = targetChain;
        mPloidyMatchType = PM_NONE;

        mShortestDistance = Links.stream().mapToLong(x -> x.length()).min().getAsLong();
    }

    public final SvChain targetChain() { return mChainTarget; }

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
            mBreakendPloidyMatched.put(breakend, Doubles.equal(bePloidy, mPloidy));
        }
    }

    public void addComDupBreakends(
            final SvBreakend compDupStart, final SvBreakend compDupEnd, double compDupPloidy,
            final SvBreakend otherStart, double otherPloidyStart,
            final SvBreakend otherEnd, double otherPloidyEnd, double otherRelativePloidy, double otherUncertainty)
    {
        if(compDupPloidy == 0 || otherPloidyStart == 0 || otherPloidyEnd== 0)
            return;

        if (copyNumbersEqual(compDupPloidy * 2, otherRelativePloidy))
        {
            mPloidyMatchType = PM_MATCHED;
        }
        else if (ploidyOverlap(compDupPloidy * 2, compDupStart.ploidyUncertainty(), otherRelativePloidy, otherUncertainty))
        {
            mPloidyMatchType = PM_OVERLAP;
        }

        if(mPloidyMatchType == PM_NONE)
            mPloidy = min(compDupPloidy, otherRelativePloidy/2);
        else
            mPloidy = (compDupPloidy + otherRelativePloidy/2) * 0.5;

        mBreakendPloidy.put(compDupStart, compDupPloidy);
        mBreakendPloidy.put(compDupEnd, compDupPloidy);
        mBreakendPloidyMatched.put(compDupStart, true);
        mBreakendPloidyMatched.put(compDupEnd, true);

        mBreakendPloidy.put(otherStart, otherPloidyStart);
        mBreakendPloidy.put(otherEnd, otherPloidyEnd);
        mBreakendPloidyMatched.put(otherStart, false);
        mBreakendPloidyMatched.put(otherEnd, false);
    }

    public void addFoldbackBreakends(
            final SvBreakend foldbackStart, final SvBreakend foldbackEnd, double foldbackPloidy,
            final SvBreakend otherBreakend, double otherPloidy, double otherRelativePloidy, double otherUncertainty)
    {
        if(foldbackPloidy == 0 || otherPloidy == 0)
            return;

        if (copyNumbersEqual(foldbackPloidy * 2, otherRelativePloidy))
        {
            mPloidyMatchType = PM_MATCHED;
        }
        else if (ploidyOverlap(foldbackPloidy * 2, foldbackStart.ploidyUncertainty(), otherRelativePloidy, otherUncertainty))
        {
            mPloidyMatchType = PM_OVERLAP;
        }

        if(mPloidyMatchType == PM_NONE)
            mPloidy = min(foldbackPloidy, otherRelativePloidy/2);
        else
            mPloidy = (foldbackPloidy + otherRelativePloidy/2) * 0.5;

        mBreakendPloidy.put(foldbackStart, foldbackPloidy);
        mBreakendPloidy.put(foldbackEnd, foldbackPloidy);
        mBreakendPloidyMatched.put(foldbackStart, true);
        mBreakendPloidyMatched.put(foldbackEnd, true);

        mBreakendPloidy.put(otherBreakend, otherPloidy);
        mBreakendPloidyMatched.put(otherBreakend, Doubles.equal(otherPloidy, mPloidy));
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

        mPloidy = min(ploidy1, ploidy2);

        if(mPloidyMatchType == PM_NONE)
        {
            mBreakendPloidyMatched.put(breakend1, Doubles.equal(ploidy1, mPloidy));
            mBreakendPloidyMatched.put(breakend2, Doubles.equal(ploidy2, mPloidy));
        }
        else
        {
            // no longer taking the average since no only is this usually above the ploidy of one of the breakends, but can
            // be above the max ploidy too
            // mPloidy = (ploidy1 + ploidy2) * 0.5;
            mBreakendPloidyMatched.put(breakend1, true);
            mBreakendPloidyMatched.put(breakend2, true);
        }
    }

    public double breakendPloidy(final SvBreakend breakend)
    {
        return mBreakendPloidy.get(breakend);
    }

    public boolean breakendPloidyMatched(final SvBreakend breakend)
    {
        return mBreakendPloidyMatched.get(breakend);
    }

    public void overrideBreakendPloidyMatched(final SvBreakend breakend)
    {
        mBreakendPloidyMatched.put(breakend, false);
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

    public ChainingRule topRule()
    {
        return mRules.get(0);
    }

    public long shortestLinkDistance() { return mShortestDistance; }

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

        if(mBreakendPloidyMatched.isEmpty() || mBreakendPloidy.isEmpty())
            return false;

        if(Links.size() == 1 && mBreakendPloidy.size() != 2)
            return false;

        if(Links.size() == 2 && mBreakendPloidy.size() < 3)
        {
            if(Links.get(0) != Links.get(1))
                return false;
        }

        return mPloidy > 0;
    }

}
