package com.hartwig.hmftools.linx.chaining;

import static com.hartwig.hmftools.linx.analysis.SvUtilities.copyNumbersEqual;
import static com.hartwig.hmftools.linx.analysis.SvUtilities.formatJcn;
import static com.hartwig.hmftools.linx.chaining.ChainJcnLimits.jcnOverlap;
import static com.hartwig.hmftools.linx.chaining.ChainingRule.COMP_DUP_SPLIT;
import static com.hartwig.hmftools.linx.chaining.ChainingRule.FOLDBACK_SPLIT;
import static com.hartwig.hmftools.linx.chaining.ChainingRule.calcRulePriority;

import java.util.List;
import java.util.Map;

import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import com.hartwig.hmftools.common.utils.Doubles;
import com.hartwig.hmftools.linx.types.LinkedPair;
import com.hartwig.hmftools.linx.types.SvBreakend;

public class ProposedLinks
{
    public final List<LinkedPair> Links;

    private final List<ChainingRule> mRules;
    private int mPriority;

    // for complex types, the chain which is being replicated by the foldback or complex dup
    private SvChain mChainTarget;

    // chain operating as a foldback - only relevant for foldback-splitting links
    private SvChain mFoldbackChain;

    private double mJcn; // the min ploidy if the breakends differ, otherwise the average
    private Map<SvBreakend, Double> mBreakendJcn;

    // indication of whether a breakend is exhausted by making these links
    private Map<SvBreakend, Boolean> mExhaustBreakend;

    private String mJcnMatchType;
    private int mShortestDistance;

    public static final String PM_MATCHED = "Matched";
    public static final String PM_NONE = "None";
    public static final String PM_OVERLAP = "Overlap";

    public ProposedLinks(LinkedPair link, ChainingRule rule)
    {
        Links = Lists.newArrayList(link);

        mRules = Lists.newArrayList(rule);
        mPriority = calcRulePriority(mRules);

        mJcn = 0;
        mBreakendJcn = Maps.newHashMap();
        mExhaustBreakend = Maps.newHashMap();
        mChainTarget = null;
        mFoldbackChain = null;
        mJcnMatchType = PM_NONE;

        mShortestDistance = link.baseLength();
    }

    public ProposedLinks(List<LinkedPair> links, ChainingRule rule, final SvChain targetChain, final SvChain foldbackChain)
    {
        Links = links;

        mRules = Lists.newArrayList(rule);
        mPriority = calcRulePriority(mRules);

        mJcn = 0;
        mBreakendJcn = Maps.newHashMap();
        mExhaustBreakend = Maps.newHashMap();
        mChainTarget = targetChain;
        mFoldbackChain = foldbackChain;
        mJcnMatchType = PM_NONE;

        mShortestDistance = Links.stream().mapToInt(x -> x.positionDistance()).min().getAsInt();
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

    public final boolean multiConnection() { return Links.size() > 1; }

    public double jcn() { return mJcn; }

    public void setLowerJcn(double jcn)
    {
        if(jcn >= mJcn || multiConnection())
            return;

        mJcn = jcn;

        // this may change whether a link is exhausted
        for(Map.Entry<SvBreakend,Double> entry : mBreakendJcn.entrySet())
        {
            final SvBreakend breakend = entry.getKey();
            double beJcn = entry.getValue();
            mExhaustBreakend.put(breakend, Doubles.equal(beJcn, mJcn));
        }
    }

    public void addComDupBreakends(
            final SvBreakend compDupStart, final SvBreakend compDupEnd, double compDupJcn,
            final SvBreakend otherStart, final SvBreakend otherEnd, double otherJcn)
    {
        if(compDupJcn == 0 || otherJcn == 0)
            return;

        // no match for this type of connection

        mJcn = compDupJcn;

        mBreakendJcn.put(compDupStart, compDupJcn);
        mBreakendJcn.put(compDupEnd, compDupJcn);
        mExhaustBreakend.put(compDupStart, true);
        mExhaustBreakend.put(compDupEnd, true);

        mBreakendJcn.put(otherStart, otherJcn);
        mBreakendJcn.put(otherEnd, otherJcn);
    }

    public void addFoldbackBreakends(
            final SvBreakend foldbackStart, final SvBreakend foldbackEnd, double foldbackJcn,
            final SvBreakend otherBreakend, double otherJcn, double otherUncertainty)
    {
        if(foldbackJcn == 0 || otherJcn == 0)
            return;

        if(copyNumbersEqual(foldbackJcn * 2, otherJcn))
        {
            mJcnMatchType = PM_MATCHED;
        }
        else if(jcnOverlap(foldbackJcn * 2, foldbackStart.jcnUncertainty(), otherJcn, otherUncertainty))
        {
            mJcnMatchType = PM_OVERLAP;
        }

        if(mJcnMatchType == PM_NONE)
            mJcn = foldbackJcn;
        else
            mJcn = (foldbackJcn + otherJcn/2) * 0.5;

        mBreakendJcn.put(foldbackStart, foldbackJcn);
        mBreakendJcn.put(foldbackEnd, foldbackJcn);

        if(foldbackStart.getSV() == foldbackEnd.getSV() && copyNumbersEqual(foldbackStart.jcn(), foldbackJcn))
        {
            // mark the foldback breakends as exhausted only if clear that they are fully used in this pair of links
            mExhaustBreakend.put(foldbackStart, true);
            mExhaustBreakend.put(foldbackEnd, true);
        }

        mBreakendJcn.put(otherBreakend, otherJcn);
    }

    public void addBreakendPloidies(
            final SvBreakend breakend1, double jcn1,
            final SvBreakend breakend2, double jcn2)
    {
        if(jcn1 == 0 || jcn2 == 0)
            return;

        mBreakendJcn.put(breakend1, jcn1);
        mBreakendJcn.put(breakend2, jcn2);

        if(copyNumbersEqual(jcn1, jcn2))
        {
            mJcnMatchType = PM_MATCHED;
        }
        else if(jcnOverlap(jcn1, breakend1.jcnUncertainty(), jcn2, breakend2.jcnUncertainty()))
        {
            mJcnMatchType = PM_OVERLAP;
        }

        if(mJcnMatchType == PM_NONE)
        {
            if(jcn1 > jcn2)
            {
                mJcn = jcn2;
            }
            else
            {
                mJcn = jcn1;
            }
        }
        else
        {
            mJcn = (jcn1 + jcn2) * 0.5;
        }
    }

    public double breakendJcn(final SvBreakend breakend)
    {
        return mBreakendJcn.get(breakend);
    }

    public boolean exhaustBreakend(final SvBreakend breakend)
    {
        Boolean exhausted = mExhaustBreakend.get(breakend);
        return exhausted != null && exhausted;
    }

    public void overrideBreakendJcnMatched(final SvBreakend breakend, boolean matched)
    {
        mExhaustBreakend.put(breakend, matched);
    }

    public final String jcnMatchType() { return mJcnMatchType; }
    public boolean linkJcnMatch() { return mJcnMatchType != PM_NONE; }

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
        return String.format("pair(%s) rule(%s) jcn(%s) match(%s) length(%d) priority(%d)",
                Links.get(0).toString(), mRules.get(0), formatJcn(mJcn),
                mJcnMatchType, mShortestDistance, mPriority);
    }

    public boolean isValid()
    {
        if(Links.isEmpty() || Links.size() > 2)
            return false;

        if(mBreakendJcn.isEmpty())
            return false;

        if(Links.size() == 1 && mBreakendJcn.size() != 2)
            return false;

        if(Links.size() == 2 && mBreakendJcn.size() < 3)
        {
            if(!Links.get(0).matches(Links.get(1)))
                return false;
        }

        return mJcn > 0;
    }

}
