package com.hartwig.hmftools.linx.chaining;

import static java.lang.Math.abs;
import static java.lang.Math.max;
import static java.lang.Math.min;

import static com.hartwig.hmftools.common.variant.structural.StructuralVariantType.INV;
import static com.hartwig.hmftools.linx.analysis.SvUtilities.copyNumbersEqual;
import static com.hartwig.hmftools.linx.analysis.SvUtilities.formatPloidy;
import static com.hartwig.hmftools.linx.chaining.ChainPloidyLimits.ploidyMatch;
import static com.hartwig.hmftools.linx.chaining.ChainPloidyLimits.ploidyOverlap;
import static com.hartwig.hmftools.linx.chaining.ChainingRule.ADJACENT;
import static com.hartwig.hmftools.linx.chaining.ChainingRule.ADJACENT_MATCH;
import static com.hartwig.hmftools.linx.chaining.ChainingRule.AP_SUPPORT;
import static com.hartwig.hmftools.linx.chaining.ChainingRule.COMP_DUP_SPLIT;
import static com.hartwig.hmftools.linx.chaining.ChainingRule.FOLDBACK;
import static com.hartwig.hmftools.linx.chaining.ChainingRule.FOLDBACK_SPLIT;
import static com.hartwig.hmftools.linx.chaining.ChainingRule.NEAREST;
import static com.hartwig.hmftools.linx.chaining.ChainingRule.PLOIDY_MATCH;
import static com.hartwig.hmftools.linx.chaining.ChainingRule.PLOIDY_MAX;
import static com.hartwig.hmftools.linx.chaining.ChainingRule.PLOIDY_OVERLAP;
import static com.hartwig.hmftools.linx.chaining.ChainingRule.ONLY;
import static com.hartwig.hmftools.linx.chaining.FoldbackBreakendPair.addByPloidy;
import static com.hartwig.hmftools.linx.chaining.FoldbackBreakendPair.containsBreakendPair;
import static com.hartwig.hmftools.linx.chaining.FoldbackBreakendPair.removeBreakendPair;
import static com.hartwig.hmftools.linx.chaining.ProposedLinks.PM_MATCHED;
import static com.hartwig.hmftools.linx.chaining.ProposedLinks.PM_NONE;
import static com.hartwig.hmftools.linx.chaining.ProposedLinks.PM_OVERLAP;
import static com.hartwig.hmftools.linx.types.SvLinkedPair.hasLinkClash;
import static com.hartwig.hmftools.linx.types.SvVarData.SE_END;
import static com.hartwig.hmftools.linx.types.SvVarData.SE_START;
import static com.hartwig.hmftools.linx.types.SvVarData.isStart;

import java.util.List;
import java.util.Map;
import java.util.stream.Collectors;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.linx.types.SvBreakend;
import com.hartwig.hmftools.linx.types.SvLinkedPair;
import com.hartwig.hmftools.linx.types.SvVarData;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;

public class ChainRuleSelector
{
    private ChainLinkAllocator mLinkAllocator;

    private int mClusterId;
    private boolean mHasReplication;
    private List<ChainingRule> mRulesToApply;

    private List<FoldbackBreakendPair> mFoldbackBreakendPairs; // a pair of breakends forming a breakend, either chained or single
    private boolean mFoldbacksInitialised;

    // references from chain-finder
    private final Map<SvBreakend, List<SvLinkedPair>> mSvBreakendPossibleLinks;
    private final ChainPloidyLimits mPloidyLimits;
    private final Map<SvVarData, SvChainState> mSvConnectionsMap;
    private final Map<SvVarData,List<SvLinkedPair>> mComplexDupCandidates;
    private final List<SvVarData> mFoldbacks;
    private final List<SvLinkedPair> mAdjacentMatchingPairs;
    private final List<SvLinkedPair> mAdjacentPairs;
    private final List<SvChain> mChains;

    private static final Logger LOGGER = LogManager.getLogger(ChainRuleSelector.class);

    public ChainRuleSelector(
            final ChainLinkAllocator linkAllocator,
                final ChainPloidyLimits ploidyLimits,
            final Map<SvBreakend, List<SvLinkedPair>> svBreakendPossibleLinks,
            final List<SvVarData> foldbacks,
            final Map<SvVarData,List<SvLinkedPair>> complexDupCandidates,
            final List<SvLinkedPair> adjacentMatchingPairs,
            final List<SvLinkedPair> adjacentPairs,
            final List<SvChain> chains)
    {
        mLinkAllocator = linkAllocator;
        mPloidyLimits = ploidyLimits;
        mSvBreakendPossibleLinks = svBreakendPossibleLinks;
        mSvConnectionsMap = linkAllocator.getSvConnectionsMap();
        mFoldbacks = foldbacks;
        mComplexDupCandidates = complexDupCandidates;
        mAdjacentMatchingPairs = adjacentMatchingPairs;
        mAdjacentPairs = adjacentPairs;
        mChains = chains;
        mRulesToApply = Lists.newArrayList();
        mFoldbackBreakendPairs = Lists.newArrayList();
        mFoldbacksInitialised = false;
    }

    public void initialise(int clusterId, boolean clusterHasReplication)
    {
        mHasReplication = clusterHasReplication;
        mClusterId = clusterId;

        mRulesToApply.clear();

        if(mHasReplication)
        {
            mRulesToApply.add(FOLDBACK_SPLIT);
            mRulesToApply.add(COMP_DUP_SPLIT);
            mRulesToApply.add(ONLY);

            mRulesToApply.add(PLOIDY_MATCH);
            mRulesToApply.add(ADJACENT);
            mRulesToApply.add(PLOIDY_MAX);
        }
        else
        {
            mRulesToApply.add(ONLY);
            mRulesToApply.add(ADJACENT_MATCH);
        }

        mRulesToApply.add(NEAREST);

        mFoldbacksInitialised = false;
        mFoldbackBreakendPairs.clear();
    }

    public List<ProposedLinks> findProposedLinks()
    {
        // find the next set of possible links to make according to the priority scheme
        // which is expressed in the set of chaining rules (ie the enumerated type)

        // some special cases
        List<ProposedLinks> proposedLinks = Lists.newArrayList();

        for(int i = 0; i < mRulesToApply.size(); ++i)
        {
            final ChainingRule rule = mRulesToApply.get(i);

            switch (rule)
            {
                case ONLY:
                    proposedLinks = findSingleOptionPairs(proposedLinks);
                    break;

                case FOLDBACK_SPLIT:
                    proposedLinks = findFoldbackLinks(proposedLinks);
                    break;

                case COMP_DUP_SPLIT:
                    proposedLinks = findComplexDupSplits(proposedLinks);
                    break;

                case PLOIDY_MATCH:
                    proposedLinks = findPloidyMatchPairs(proposedLinks);
                    break;

                case ADJACENT:
                    proposedLinks = findAdjacentPairs(proposedLinks);
                    break;

                case ADJACENT_MATCH:
                    proposedLinks = findAdjacentMatchingPairs(proposedLinks);
                    break;

                case PLOIDY_MAX:
                    proposedLinks = findHighestPloidy(proposedLinks);
                    break;

                case NEAREST:
                    proposedLinks = findNearest(proposedLinks);
                    break;
            }

            mLinkAllocator.removeSkippedPairs(proposedLinks);

            // the top-priority rule acts like an annotation - whether the link violates the cluster ploidy across
            // one or more breakend segments - so even if a higher rule has only 1 proposed link, if it violates this condition
            // then keep searching for one which doesn't
            boolean linksHavePloidySupport = checkClusterPloidySupport(proposedLinks);

            // the last rule 'nearest' can throw up multiple possible links but since they're not conflicting and ordered
            // from shortest to longest, there's no need to cull any more
            if(i > 0 && rule != NEAREST)
                cullByPriority(proposedLinks);

            // if the proposed links have been reduced to a single-top priority rule, which also meets the ploidy-support restriction
            // then take it
            if(proposedLinks.size() == 1 && linksHavePloidySupport)
                return proposedLinks;
        }

        return proposedLinks;
    }

    private List<ProposedLinks> findSingleOptionPairs(List<ProposedLinks> proposedLinks)
    {
        // find all breakends with only one other link option
        if(!proposedLinks.isEmpty())
        {
            for(final ProposedLinks proposedLink : proposedLinks)
            {
                boolean hasSingleOption = false;

                for(final SvLinkedPair pair : proposedLink.Links)
                {
                    for(int se = SE_START; se <= SE_END; ++se)
                    {
                        List<SvLinkedPair> possiblePairs = mSvBreakendPossibleLinks.get(pair.getBreakend(isStart(se)));

                        if(possiblePairs == null)
                        {
                            LOGGER.error("breakend({}) has no possible pairs, from proposedLink: {}",
                                    pair.getBreakend(isStart(se)), proposedLink);
                            continue;
                        }

                        if(possiblePairs.size() == 1)
                        {
                            hasSingleOption = true;
                            break;
                        }
                    }

                    if(hasSingleOption)
                    {
                        proposedLink.addRule(ONLY);
                        break;
                    }
                }
            }
        }

        for(Map.Entry<SvBreakend, List<SvLinkedPair>> entry : mSvBreakendPossibleLinks.entrySet())
        {
            if(entry.getValue().isEmpty())
            {
                LOGGER.warn("breakend({}) has no possibles left, should be purged", entry.getKey().toString());
                continue;
            }

            if(entry.getValue().size() >= 2) // disable connections to an INV for now
                continue;

            SvBreakend limitingBreakend = entry.getKey();

            // isSpecificSV(limitingBreakend.getSV());

            final SvLinkedPair newPair;

            if(entry.getValue().size() == 2)
            {
                // consider a link only to an INV as a single option
                final SvLinkedPair pair1 = entry.getValue().get(0);
                final SvLinkedPair pair2 = entry.getValue().get(1);
                final SvVarData otherSv1 = pair1.getOtherSV(limitingBreakend.getSV());
                final SvVarData otherSv2 = pair2.getOtherSV(limitingBreakend.getSV());

                if(otherSv1 == otherSv2 && otherSv1.type() == INV)
                {
                    newPair = pair1.length() < pair2.length() ? pair1 : pair2;
                }
                else
                {
                    continue;
                }
            }
            else
            {
                newPair = entry.getValue().get(0);
            }

            if(mLinkAllocator.hasSkippedPairs(newPair))
                continue;

            // skip the duplicate link stored against the other breakend
            if(proposedLinks.stream().map(x -> x.Links.get(0)).anyMatch(y -> y == newPair))
                continue;

            double ploidyFirst = mLinkAllocator.getUnlinkedBreakendCount(newPair.firstBreakend(), true);
            double ploidySecond = mLinkAllocator.getUnlinkedBreakendCount(newPair.secondBreakend(), true);

            if(ploidyFirst == 0 || ploidySecond == 0)
                continue;

            ProposedLinks proposedLink = new ProposedLinks(newPair, ONLY);
            proposedLink.addBreakendPloidies(newPair.firstBreakend(), ploidyFirst, newPair.secondBreakend(), ploidySecond);

            // check for another proposed link with a clashing breakend, and if found take the lower ploidy and short link
            boolean addNew = true;
            int index = 0;
            while(index < proposedLinks.size())
            {
                final ProposedLinks otherLink = proposedLinks.get(index);
                final SvLinkedPair otherPair = otherLink.Links.get(0);

                if(!otherPair.hasLinkClash(newPair) && !otherPair.oppositeMatch(newPair))
                {
                    ++index;
                    continue;
                }

                if(mHasReplication)
                {
                    if (copyNumbersEqual(otherLink.ploidy(), proposedLink.ploidy()))
                    {
                        if (proposedLink.shortestLinkDistance() > otherLink.shortestLinkDistance())
                        {
                            addNew = false;
                            break;
                        }

                        proposedLinks.remove(otherLink);
                    }

                    // keep both for now and let downstream rules decide
                    ++index;
                }
                else
                {
                    if (proposedLink.shortestLinkDistance() > otherLink.shortestLinkDistance())
                    {
                        addNew = false;
                        break;
                    }

                    proposedLinks.remove(otherLink);

                }
            }

            if(addNew)
                proposedLinks.add(proposedLink);
        }

        return proposedLinks;
    }

    private void updateFoldbackBreakends()
    {
        List<FoldbackBreakendPair> existingChainedPairs = Lists.newArrayList();

        if(!mFoldbacksInitialised)
        {
            // add non-chained simple INV foldbacks first - chained ones will be added next
            mFoldbacksInitialised = true;
            mFoldbackBreakendPairs.clear();

            for (SvVarData foldback : mFoldbacks)
            {
                if (foldback.isSingleBreakendFoldback() || foldback.isChainedFoldback())
                    continue;

                SvBreakend foldbackStart = foldback.getBreakend(true);
                SvBreakend foldbackEnd = foldback.getBreakend(false);


                double ploidy = min(mLinkAllocator.getUnlinkedBreakendCount(foldbackStart), mLinkAllocator.getUnlinkedBreakendCount(foldbackEnd));

                if (ploidy == 0)
                    continue;

                FoldbackBreakendPair fbPair = new FoldbackBreakendPair(foldbackStart, foldbackEnd, ploidy, null);

                addByPloidy(mFoldbackBreakendPairs, fbPair);
            }
        }
        else
        {
            // remove any foldbacks with exhausted breakends
            int index = 0;
            while(index < mFoldbackBreakendPairs.size())
            {
                final FoldbackBreakendPair fbPair = mFoldbackBreakendPairs.get(index);

                double minUnlinkedPloidy = min(
                        mLinkAllocator.getUnlinkedBreakendCount(fbPair.BreakendStart, true),
                        mLinkAllocator.getUnlinkedBreakendCount(fbPair.BreakendEnd, true));

                if(minUnlinkedPloidy == 0)
                {
                    LOGGER.debug("foldback pair({}) removed from consideration", fbPair);
                    mFoldbackBreakendPairs.remove(index);
                    continue;
                }

                if(fbPair.isChained())
                    existingChainedPairs.add(fbPair); // will be checked against current chains next up
                else
                    fbPair.Ploidy = minUnlinkedPloidy;

                ++index;
            }
        }

        // also include chain ends which effectively form a foldback
        for(final SvChain chain : mChains)
        {
            final SvBreakend chainStart = chain.getOpenBreakend(true);
            final SvBreakend chainEnd = chain.getOpenBreakend(false);

            if(chainStart == null || chainEnd == null)
                continue;

            // are breakends consecutive and facing the same direction?
            if(!chainStart.chromosome().equals(chainEnd.chromosome()) || chainStart.orientation() != chainEnd.orientation())
                continue;

            int startIndex = chainStart.getClusterChrPosIndex();
            int endIndex = chainEnd.getClusterChrPosIndex();

            if(chainStart != chainEnd && abs(startIndex - endIndex) != 1)
            {
                if(abs(startIndex - endIndex) > 2)
                    continue;

                // the foldback is invalid if it has a deletion bridge with overhang on the front-facing breakend
                boolean startIsFront = (chainStart.position() < chainEnd.position()) == (chainStart.orientation() == 1);
                final SvBreakend frontBreakend = startIsFront ? chainStart : chainEnd;

                // check for a DB which either allows or invalidates this foldback
                final SvLinkedPair dbLink = frontBreakend.getSV().getDBLink(frontBreakend.usesStart());

                if(dbLink != null && dbLink.length() > 0)
                    continue;
            }

            FoldbackBreakendPair fbPair = new FoldbackBreakendPair(chainStart, chainEnd, chain.ploidy(), chain);

            removeBreakendPair(existingChainedPairs, fbPair);

            if(!containsBreakendPair(mFoldbackBreakendPairs, fbPair))
            {
                LOGGER.debug("chain({}) adding chained foldback breakends({})", chain.id(), fbPair);
                addByPloidy(mFoldbackBreakendPairs, fbPair);
            }
        }

        // remove any chained foldbacks no longer supported by the current set of chains
        existingChainedPairs.stream().forEach(x -> removeBreakendPair(mFoldbackBreakendPairs, x));
    }

    private List<ProposedLinks> findFoldbackLinks(final List<ProposedLinks> proposedLinks)
    {
        updateFoldbackBreakends();

        if (mFoldbackBreakendPairs.isEmpty())
            return proposedLinks;

        // find the highest ploidy foldback where:
        //	a) it splits another chain with 2x ploidy
        //	b) it matches the ploidy of another chain
        //	c) it splits another foldback with >2x ploidy
        //	d) it is itself split by another foldback or complex duplication

        List<ProposedLinks> newProposedLinks = Lists.newArrayList();
        int linkScore = -1; // keep track of the highest priority proposed links amongst the foldbacks

        double lastFoldbackPloidy = 0;

        // foldbacks are cached from highest to lowest ploidy already
        for (final FoldbackBreakendPair fbPair : mFoldbackBreakendPairs)
        {
            // break if the next foldback pair has a lower ploidy than any of those which have proposed a link
            if (!newProposedLinks.isEmpty() && !copyNumbersEqual(fbPair.Ploidy, lastFoldbackPloidy))
                break;

            lastFoldbackPloidy = fbPair.Ploidy;
            SvBreakend foldbackStart = fbPair.BreakendStart;
            SvBreakend foldbackEnd = fbPair.BreakendEnd;

            SvVarData foldback = foldbackStart.getSV();

            // if a foldback breakend is chained, use the ploidy of its chain
            SvChain foldbackChain = null;
            double foldbackPloidy = fbPair.Ploidy;
            double foldbackUncertainty = 0;

            if (fbPair.isChained())
            {
                foldbackChain = fbPair.Chain;
                foldbackUncertainty = foldbackChain.ploidyUncertainty();
            }
            else
            {
                foldback.ploidyUncertainty();
            }

            // search through all available pairs for the top priority type of link for this foldback
            List<SvLinkedPair> pairsOnFbStart = mSvBreakendPossibleLinks.get(foldbackStart);
            List<SvLinkedPair> pairsOnFbEnd = mSvBreakendPossibleLinks.get(foldbackEnd);

            if (pairsOnFbStart == null || pairsOnFbEnd == null)
                continue;

            pairsOnFbStart = Lists.newArrayList(pairsOnFbStart);
            pairsOnFbEnd = Lists.newArrayList(pairsOnFbEnd );

            cullDualOptionPairs(foldbackStart, pairsOnFbStart);

            for (SvLinkedPair pairStart : pairsOnFbStart)
            {
                SvVarData nonFbVar = pairStart.getOtherSV(foldback);
                SvBreakend otherBreakend = pairStart.getOtherBreakend(foldbackStart);

                // find the other pairing - would expect this to exist
                SvLinkedPair pairEnd = pairsOnFbEnd.stream()
                        .filter(x -> x.getOtherBreakend(foldbackEnd) == otherBreakend)
                        .findFirst().orElse(null);

                if (pairEnd == null)
                    continue;

                // first establish the available ploidy of this breakend and whether it's chained
                final BreakendPloidy nonFbPloidyData = mLinkAllocator.getBreakendPloidyData(otherBreakend);
                double nonFbPloidy = nonFbPloidyData.unlinkedPloidy();

                // for low ploidy SVs, ploidy comparisons are imprecise, so don't allow splits
                boolean allowSplits = !(copyNumbersEqual(foldbackPloidy, nonFbPloidy) && foldbackPloidy < 1 && nonFbPloidy < 1);

                // a) first check if it splits another chain or SV with 2x ploidy
                if (allowSplits && ploidyMatch(foldbackPloidy * 2, foldbackUncertainty, nonFbPloidy, nonFbVar.ploidyUncertainty())
                && !nonFbPloidyData.multiConnections())
                {
                    // a 2:1 splitting event
                    if (linkScore < 3)
                    {
                        newProposedLinks.clear();
                        linkScore = 3;
                    }

                    ProposedLinks proposedLink = new ProposedLinks(
                            Lists.newArrayList(pairStart, pairEnd), FOLDBACK_SPLIT, nonFbPloidyData.MaxPloidyChain, foldbackChain);

                    proposedLink.addFoldbackBreakends(
                            foldbackStart, foldbackEnd, foldbackPloidy,
                            otherBreakend, nonFbPloidy, nonFbVar.ploidyUncertainty());

                    LOGGER.trace("type-A: foldback breakends({} & {}) ploidy({}) exact split of breakend({}) ploidy({})",
                            foldbackStart, foldbackEnd, formatPloidy(foldbackPloidy), otherBreakend, formatPloidy(nonFbPloidy));

                    newProposedLinks.add(proposedLink);
                    continue;
                }

                if (linkScore > 2)
                    continue;

                // b) check for an exact match with a chain or another SV
                if (!nonFbPloidyData.multiConnections() && ploidyMatch(foldbackPloidy, foldbackUncertainty, nonFbPloidy, nonFbVar.ploidyUncertainty()))
                {
                    if (linkScore < 3)
                    {
                        newProposedLinks.clear();
                        linkScore = 3;
                    }

                    // where there is a choice to be made between which breakends are used as is usually the case of foldbacks,
                    // choose the breakend with the most unlinked ploidy, which implies its other end is already chained
                    final SvChainState fbConn = mLinkAllocator.getSvConnectionsMap().get(foldback);

                    final SvBreakend fbBreakend = fbConn.unlinked(foldbackStart.usesStart()) > fbConn.unlinked(foldbackEnd.usesStart())
                            ? foldbackStart : foldbackEnd;

                    ProposedLinks proposedLink = new ProposedLinks(SvLinkedPair.from(fbBreakend, otherBreakend), FOLDBACK);
                    proposedLink.addBreakendPloidies(fbBreakend, foldbackPloidy, otherBreakend, nonFbPloidy);

                    LOGGER.trace("type-B: foldback({}) ploidy({}) matched with breakend({}) ploidy({})",
                            fbBreakend, formatPloidy(foldbackPloidy), otherBreakend, formatPloidy(nonFbPloidy));

                    newProposedLinks.add(proposedLink);
                    continue;
                }

                if (linkScore > 1)
                    continue;

                // c) check whether this foldback splits another foldback with > 2x ploidy
                FoldbackBreakendPair otherFbPair = mFoldbackBreakendPairs.stream()
                        .filter(x -> x != fbPair)
                        .filter(x -> x.BreakendStart == otherBreakend || x.BreakendEnd == otherBreakend)
                        .findFirst().orElse(null);

                // if this foldback splits another foldback then the other foldback's ploidy must be 2x or higher
                if (otherFbPair != null
                        && (!ploidyMatch(foldbackPloidy * 2, foldbackUncertainty, nonFbPloidy, nonFbVar.ploidyUncertainty())
                        && nonFbPloidy > foldbackPloidy * 2))
                {
                    if (linkScore < 1)
                    {
                        newProposedLinks.clear();
                        linkScore = 1;
                    }

                    ProposedLinks proposedLink = new ProposedLinks(
                            Lists.newArrayList(pairStart, pairEnd), FOLDBACK_SPLIT, nonFbPloidyData.MaxPloidyChain, foldbackChain);

                    proposedLink.addFoldbackBreakends(
                            foldbackStart, foldbackEnd, foldbackPloidy,
                            otherBreakend, nonFbPloidy, nonFbVar.ploidyUncertainty());

                    newProposedLinks.add(proposedLink);

                    LOGGER.trace("type-C: foldback breakends({} & {}) ploidy({}) non-exact split of foldback breakend({}) ploidy({})",
                            foldbackStart, foldbackEnd, formatPloidy(foldbackPloidy), otherBreakend, formatPloidy(nonFbPloidy));

                    continue;
                }

                // d) check whether the foldback is itself split by another foldback or complex duplication
                if (otherFbPair != null && allowSplits
                && ploidyMatch(foldbackPloidy, foldbackUncertainty, nonFbPloidy * 2, nonFbVar.ploidyUncertainty()))
                {
                    linkScore = 0;

                    // a 2:1 splitting event in reverse
                    ProposedLinks proposedLink = new ProposedLinks(
                            Lists.newArrayList(
                                    SvLinkedPair.from(otherFbPair.BreakendStart, foldbackStart),
                                    SvLinkedPair.from(otherFbPair.BreakendEnd, foldbackStart)),
                            FOLDBACK_SPLIT, foldbackChain, otherFbPair.Chain);

                    proposedLink.addFoldbackBreakends(
                            otherFbPair.BreakendStart, otherFbPair.BreakendEnd, otherFbPair.Ploidy,
                            foldbackStart, foldbackPloidy, foldbackUncertainty);

                    newProposedLinks.add(proposedLink);

                    LOGGER.trace("type-D: foldback breakend({}) ploidy({}) split by other foldback pair({})",
                            foldbackStart, formatPloidy(foldbackPloidy), otherFbPair);

                    continue;
                }
            }
        }

            return newProposedLinks;
    }

    private void cullDualOptionPairs(final SvBreakend sourceBreakend, List<SvLinkedPair> pairs)
    {
        // if this breakend can connect to both ends of another SV (ie an INV), then removed the end which
        // is most exhausted (that is most chained)
        int i = 0;
        while(i < pairs.size())
        {
            final SvBreakend breakend = pairs.get(i).getOtherBreakend(sourceBreakend);

            boolean breakendRemoved = false;
            for(int j = i + 1; j < pairs.size(); ++j)
            {
                final SvBreakend breakend2 = pairs.get(j).getOtherBreakend(sourceBreakend);

                if(breakend.getSV() != breakend2.getSV())
                    continue;

                final SvChainState svConn = mLinkAllocator.getSvConnectionsMap().get(breakend.getSV());

                if(svConn == null)
                    continue;

                if(copyNumbersEqual(svConn.unlinked(breakend.usesStart()), svConn.unlinked(breakend2.usesStart()))
                && pairs.get(i).length() > pairs.get(j).length())
                {
                    pairs.remove(i);
                    breakendRemoved = true;
                }
                else if(svConn.unlinked(breakend.usesStart()) < svConn.unlinked(breakend2.usesStart()))
                {
                    pairs.remove(i);
                    breakendRemoved = true;
                }
                else
                {
                    pairs.remove(j);
                }

                break;
            }

            if(!breakendRemoved)
                ++i;
        }
    }

    private List<ProposedLinks> findComplexDupSplits(List<ProposedLinks> proposedLinks)
    {
        // both ends of a foldback or complex DUP connect to one end of another SV with ploidy >= 2x
        if(mComplexDupCandidates.isEmpty())
            return proposedLinks;

        // only search for complex dups when no higher priority links have been found
        if(proposedLinks.stream().anyMatch(x -> x.hasRule(FOLDBACK_SPLIT)))
            return proposedLinks;

        List<ProposedLinks> newProposedLinks = Lists.newArrayList();

        // the complex DUP SVs need to connect to both ends of either a single SV or a set of chained SVs with twice the ploidy
        // for a complex DUP of the form D - A - D, where A has double ploidy of D, and both ends of D connect to both ends of A
        List<SvVarData> invalidCompDups = Lists.newArrayList();

        for(Map.Entry<SvVarData,List<SvLinkedPair>> entry : mComplexDupCandidates.entrySet())
        {
            SvVarData compDup = entry.getKey();
            double compDupPloidy = mLinkAllocator.getUnlinkedCount(compDup);

            if(compDupPloidy == 0)
                continue;

            SvBreakend compDupBeStart = compDup.getBreakend(true);
            SvBreakend compDupBeEnd = compDup.getBreakend(false);

            List<SvLinkedPair> compDupPairs = entry.getValue();

            if(compDupPairs.size() != 2)
            {
                invalidCompDups.add(compDup);
                continue;
            }

            final SvLinkedPair pair1 = compDupPairs.get(0);
            final SvLinkedPair pair2 = compDupPairs.get(1);

            final SvVarData otherSv1 = pair1.getOtherSV(compDup);
            final SvVarData otherSv2 = pair2.getOtherSV(compDup);

            if (otherSv1 == otherSv2)
            {
                // a single SV duplicated by the complex DUP
                final SvBreakend otherBreakend1 = otherSv1.getBreakend(true);
                final SvBreakend otherBreakend2 = otherSv1.getBreakend(false);

                if(otherBreakend1 == null || otherBreakend2 == null)
                {
                    LOGGER.error("comp dup other breakends null");
                    continue;
                }

                double firstPloidy = mLinkAllocator.getUnlinkedBreakendCount(otherBreakend1);
                double secondPloidy = mLinkAllocator.getUnlinkedBreakendCount(otherBreakend2);

                if (firstPloidy == 0 || secondPloidy == 0
                || !copyNumbersEqual(compDupPloidy * 2, firstPloidy) && firstPloidy <= compDupPloidy)
                {
                    invalidCompDups.add(compDup);
                    continue;
                }

                ProposedLinks proposedLink = new ProposedLinks(Lists.newArrayList(pair1, pair2), COMP_DUP_SPLIT, null, null);

                proposedLink.addComDupBreakends(
                        compDupBeStart, compDupBeEnd, compDupPloidy,
                        otherBreakend1, firstPloidy, otherBreakend2, secondPloidy,
                        firstPloidy, otherBreakend2.ploidyUncertainty());

                newProposedLinks.add(proposedLink);

                LOGGER.trace("comDup({}) ploidy({}) matched with breakends({} & {}) ploidy({} & {})",
                        compDup.id(), formatPloidy(compDupPloidy), otherBreakend1,
                        otherBreakend2, formatPloidy(firstPloidy), formatPloidy(secondPloidy));
            }
            else
            {
                final SvBreakend otherBreakend1 = pair1.hasBreakend(compDupBeStart) ?
                        pair1.getOtherBreakend(compDupBeStart) : pair1.getOtherBreakend(compDupBeEnd);

                final SvBreakend otherBreakend2 = pair2.hasBreakend(compDupBeStart) ?
                        pair2.getOtherBreakend(compDupBeStart) : pair2.getOtherBreakend(compDupBeEnd);

                // search existing chains for open chain ends match the set of possibles for the complex DUP and with twice the ploidy
                for (SvChain chain : mChains)
                {
                    if (!copyNumbersEqual(compDupPloidy * 2, chain.ploidy()) && chain.ploidy() <= compDupPloidy)
                        continue;

                    SvBreakend chainBeStart = chain.getOpenBreakend(true);
                    SvBreakend chainBeEnd = chain.getOpenBreakend(false);

                    if(!((chainBeStart == otherBreakend1 && chainBeEnd == otherBreakend2)
                    || (chainBeStart == otherBreakend2 && chainBeEnd == otherBreakend1)))
                        continue;

                    double chainStartPloidy = mLinkAllocator.getUnlinkedBreakendCount(chainBeStart);
                    double chainEndPloidy = mLinkAllocator.getUnlinkedBreakendCount(chainBeEnd);

                    if (chainStartPloidy == 0 || chainEndPloidy == 0)
                        continue;

                    ProposedLinks proposedLink = new ProposedLinks(Lists.newArrayList(pair1, pair2), COMP_DUP_SPLIT, chain, null);

                    proposedLink.addComDupBreakends(
                            compDupBeStart, compDupBeEnd, compDupPloidy,
                            chainBeStart, chainStartPloidy, chainBeEnd, chainEndPloidy, chain.ploidy(), chain.ploidyUncertainty());

                    newProposedLinks.add(proposedLink);

                    LOGGER.trace("comDup({}) ploidy({}) matched with chain breakends({} & {}) ploidy({})",
                            compDup.id(), formatPloidy(compDupPloidy),
                            chainBeStart.toString(), chainBeEnd.toString(), formatPloidy(chain.ploidy()));

                    break;
                }
            }

            /*
            mDiagnostics.logCsv("COMP_DUP", compDup,
                    String.format("ploidy(%.1f-%.1f-%.1f) beStart(%s ploidy=%.1f-%.1f) beEnd(%s ploidy=%.1f-%.1f)",
                            compDup.ploidyMin(), compDup.ploidy(), compDup.ploidyMax(),
                            chainBeStart.toString(), chainBeStart.getSV().ploidyMin(), chainBeStart.getSV().ploidyMax(),
                            chainBeEnd.toString(), chainBeEnd.getSV().ploidyMin(), chainBeEnd.getSV().ploidyMax()));
            */
        }

        for(SvVarData compDup : invalidCompDups)
        {
            mComplexDupCandidates.remove(compDup);
        }

        // now check for a match between any previously identified proposed links and this set
        /*
        if(proposedLinks.isEmpty())
            return newProposedLinks;

        if(newProposedLinks.isEmpty())
            return proposedLinks;

        // only retain those foldback proposed links which are in the only-option set
        List<ProposedLinks> mergedProposed = newProposedLinks.stream()
                .filter(x -> proposedLinks.stream().anyMatch(y -> anyLinkMatch(x.Links, y.Links)))
                .collect(Collectors.toList());

        if(mergedProposed.isEmpty())
            return proposedLinks;

        mergedProposed.stream().forEach(x -> x.addRule(ONLY));

        return mergedProposed;
        */

        return restrictProposedLinks(newProposedLinks, proposedLinks, ONLY);
    }

    private List<ProposedLinks> findPloidyMatchPairs(List<ProposedLinks> proposedLinks)
    {
        // find pairs of matching ploidy breakends, taking the shortest where multiple exist
        List<ProposedLinks> newProposedLinks = Lists.newArrayList();

        boolean hasPloidySupportLinks = anyLinksHavePloidySupport(proposedLinks);

        if(!proposedLinks.isEmpty())
        {
            proposedLinks.stream().filter(x -> x.ploidyMatchType() == PM_MATCHED).forEach(x -> x.addRule(PLOIDY_MATCH));
            proposedLinks.stream().filter(x -> x.ploidyMatchType() == PM_OVERLAP).forEach(x -> x.addRule(PLOIDY_OVERLAP));

            if(hasPloidySupportLinks)
                return proposedLinks;
        }

        // double currentMaxPloidy = 0;
        List<SvLinkedPair> addedLinks = Lists.newArrayList();

        for(SvChainState svConn : mSvConnectionsMap.values())
        {
            SvVarData var = svConn.SV;

            // check whether this SV has any possible links with SVs of the same (remaining) rep count
            for(int be = SE_START; be <= SE_END; ++be)
            {
                if(var.isSglBreakend() && be == SE_END)
                    continue;

                boolean isStart = isStart(be);
                double breakendPloidy = svConn.unlinked(be);

                if(breakendPloidy == 0)
                    continue;

                // if(!copyNumbersEqual(breakendPloidy, currentMaxPloidy) && breakendPloidy < currentMaxPloidy)
                //    continue;

                final SvBreakend breakend = var.getBreakend(isStart);
                final List<SvLinkedPair> svLinks = mSvBreakendPossibleLinks.get(breakend);

                if(svLinks == null)
                    continue;

                ProposedLinks bestProposedLink = null;

                for(final SvLinkedPair pair : svLinks)
                {
                    if(mLinkAllocator.hasSkippedPairs(pair))
                        continue;

                    if(addedLinks.contains(pair))
                        continue;

                    SvBreakend otherBreakend = pair.getOtherBreakend(breakend);

                    double otherBreakendPloidy = mLinkAllocator.getUnlinkedBreakendCount(otherBreakend, true);

                    if(otherBreakendPloidy == 0)
                        continue;

                    String ploidyMatch = PM_NONE;
                    if(copyNumbersEqual(otherBreakendPloidy, breakendPloidy))
                    {
                        ploidyMatch = PM_MATCHED;
                    }
                    else if(ploidyOverlap(var.ploidyUncertainty(), breakendPloidy, otherBreakendPloidy, otherBreakend.ploidyUncertainty()))
                    {
                        ploidyMatch = PM_OVERLAP;
                    }
                    else
                    {
                        continue;
                    }

                    // LOGGER.trace("pair({}) with {} ploidy({} & {})",
                    //        pair.toString(), ploidyMatch, formatPloidy(breakendPloidy), formatPloidy(otherBreakendPloidy));

                    if(bestProposedLink != null && bestProposedLink.topRule() == PLOIDY_OVERLAP && ploidyMatch != PM_MATCHED)
                    {
                        // no better and longer so keep searching for a better match
                        continue;
                    }

                    ProposedLinks proposedLink = (ploidyMatch == PM_MATCHED) ?
                            new ProposedLinks(pair, PLOIDY_MATCH) : new ProposedLinks(pair, PLOIDY_OVERLAP);

                    proposedLink.addBreakendPloidies(breakend, breakendPloidy, otherBreakend, otherBreakendPloidy);

                    bestProposedLink = proposedLink;
                    addedLinks.add(pair);

                    if(ploidyMatch == PM_MATCHED) // no need to keep searching through this breakend
                        break;

                    // if(!copyNumbersEqual(proposedLink.ploidy(), currentMaxPloidy))
                    //    currentMaxPloidy = proposedLink.ploidy();
                }

                if(bestProposedLink != null)
                    newProposedLinks.add(bestProposedLink);
            }
        }

        checkClusterPloidySupport(newProposedLinks);

        // if a new link has ploidy support it will top anything found already (see earlier exit condition for proposed links)
        if(proposedLinks.isEmpty() || anyLinksHavePloidySupport(newProposedLinks))
            return newProposedLinks;
        else
            return proposedLinks;
    }

    private List<ProposedLinks> findAdjacentMatchingPairs(List<ProposedLinks> proposedLinks)
    {
        if(mAdjacentMatchingPairs.isEmpty())
            return proposedLinks;

        if(!proposedLinks.isEmpty())
        {
            for(ProposedLinks proposedLink : proposedLinks)
            {
                // check for any links which are in the adjacent set
                if(proposedLink.Links.stream().anyMatch(x -> mAdjacentMatchingPairs.contains(x)))
                {
                    proposedLink.addRule(ADJACENT);
                }
            }

            return proposedLinks;
        }

        List<ProposedLinks> newProposedLinks = Lists.newArrayList();

        while(!mAdjacentMatchingPairs.isEmpty())
        {
            SvLinkedPair nextPair = mAdjacentMatchingPairs.get(0);

            mAdjacentMatchingPairs.remove(0);

            if(mLinkAllocator.matchesExistingPair(nextPair))
                continue;

            double ploidyFirst = mLinkAllocator.getUnlinkedBreakendCount(nextPair.firstBreakend(), true);
            double ploidySecond = mLinkAllocator.getUnlinkedBreakendCount(nextPair.secondBreakend(), true);

            if(ploidyFirst == 0 || ploidySecond == 0)
                continue;

            if(!copyNumbersEqual(ploidyFirst, ploidySecond))
                continue;

            // take the average ploidy or calculate a weighted ploidy already?
            // if these links have already been partially used, then incorrect to calculate a weighted ploidy
            ProposedLinks proposedLink = new ProposedLinks(nextPair, ADJACENT);
            proposedLink.addRule(PLOIDY_MATCH);
            proposedLink.addBreakendPloidies(nextPair.firstBreakend(), ploidyFirst, nextPair.secondBreakend(), ploidySecond);
            newProposedLinks.add(proposedLink);
        }

        return restrictProposedLinks(proposedLinks, newProposedLinks, ADJACENT);
    }

    private List<ProposedLinks> findAdjacentPairs(List<ProposedLinks> proposedLinks)
    {
        if(mAdjacentPairs.isEmpty())
            return proposedLinks;

        boolean hasPloidySupportLinks = anyLinksHavePloidySupport(proposedLinks);

        if(!proposedLinks.isEmpty())
        {
            // check for any links which are in the adjacent set
            proposedLinks.stream()
                    .filter(x -> x.Links.stream().anyMatch(y -> mAdjacentPairs.contains(y)))
                    .forEach(x -> x.addRule(ADJACENT));

            if(hasPloidySupportLinks)
                return proposedLinks;
        }

        List<ProposedLinks> newProposedLinks = Lists.newArrayList();

        while(!mAdjacentPairs.isEmpty())
        {
            SvLinkedPair nextPair = mAdjacentPairs.get(0);

            mAdjacentPairs.remove(0);

            if(mLinkAllocator.matchesExistingPair(nextPair))
                continue;

            double ploidyFirst = mLinkAllocator.getUnlinkedBreakendCount(nextPair.firstBreakend(), true);
            double ploidySecond = mLinkAllocator.getUnlinkedBreakendCount(nextPair.secondBreakend(), true);

            if(ploidyFirst == 0 || ploidySecond == 0)
                continue;

            // take the average ploidy or calculate a weighted ploidy already?
            // if these links have already been partially used, then incorrect to calculate a weighted ploidy

            ProposedLinks proposedLink = new ProposedLinks(nextPair, ADJACENT);
            proposedLink.addRule(PLOIDY_MATCH);
            proposedLink.addBreakendPloidies(nextPair.firstBreakend(), ploidyFirst, nextPair.secondBreakend(), ploidySecond);
            newProposedLinks.add(proposedLink);
        }

        checkClusterPloidySupport(newProposedLinks);

        // if a new link has ploidy support it will top anything found already (see earlier exit condition for proposed links)
        if(proposedLinks.isEmpty() || anyLinksHavePloidySupport(newProposedLinks))
            return newProposedLinks;
        else
            return proposedLinks;

        // return restrictProposedLinks(proposedLinks, newProposedLinks, ADJACENT);
    }

    private List<ProposedLinks> findHighestPloidy(List<ProposedLinks> proposedLinks)
    {
        List<ProposedLinks> newProposedLinks = Lists.newArrayList();

        boolean hasPloidySupportLinks = anyLinksHavePloidySupport(proposedLinks);

        if(!proposedLinks.isEmpty())
        {
            // take the highest from amongst the proposed links
            double maxPloidy = proposedLinks.stream().mapToDouble(x -> x.ploidy()).max().getAsDouble();

            proposedLinks.stream().filter(x -> copyNumbersEqual(maxPloidy, x.ploidy())).forEach(x -> x.addRule(PLOIDY_MAX));

            if(hasPloidySupportLinks)
                return proposedLinks;
        }

        double currentMaxPloidy = 0;
        List<SvLinkedPair> addedLinks = Lists.newArrayList();

        for(SvChainState svConn : mSvConnectionsMap.values())
        {
            SvVarData var = svConn.SV;

            // check whether this SV has any possible links with SVs of the same (remaining) rep count
            for(int be = SE_START; be <= SE_END; ++be)
            {
                if(var.isSglBreakend() && be == SE_END)
                    continue;

                boolean isStart = isStart(be);
                double breakendPloidy = svConn.unlinked(be);

                if(breakendPloidy == 0)
                    continue;

                if(!copyNumbersEqual(breakendPloidy, currentMaxPloidy) && breakendPloidy < currentMaxPloidy)
                    continue;

                final SvBreakend breakend = var.getBreakend(isStart);
                final List<SvLinkedPair> svLinks = mSvBreakendPossibleLinks.get(breakend);

                breakendPloidy = mLinkAllocator.getUnlinkedBreakendCount(breakend, true);

                if(svLinks == null)
                    continue;

                for(final SvLinkedPair pair : svLinks)
                {
                    if(addedLinks.contains(pair))
                        continue;

                    if(mLinkAllocator.hasSkippedPairs(pair))
                        continue;

                    SvBreakend otherBreakend = pair.getOtherBreakend(breakend);

                    double otherBreakendPloidy = mLinkAllocator.getUnlinkedBreakendCount(otherBreakend, true);

                    if(otherBreakendPloidy == 0)
                        continue;

                    double minPairPloidy = min(otherBreakendPloidy, breakendPloidy);

                    if(!copyNumbersEqual(currentMaxPloidy, minPairPloidy) && minPairPloidy < currentMaxPloidy)
                        continue;

                    currentMaxPloidy = max(minPairPloidy, currentMaxPloidy);

                    LOGGER.trace("pair({}) with max ploidy({} & {})",
                            pair.toString(), formatPloidy(breakendPloidy), formatPloidy(otherBreakendPloidy));

                    ProposedLinks proposedLink = new ProposedLinks(pair, PLOIDY_MAX);
                    proposedLink.addBreakendPloidies(breakend, breakendPloidy, otherBreakend, otherBreakendPloidy);
                    newProposedLinks.add(proposedLink);
                }
            }
        }

        checkClusterPloidySupport(newProposedLinks);

        // if a new link has ploidy support it will top anything found already (see earlier exit condition for proposed links)
        if(proposedLinks.isEmpty() || anyLinksHavePloidySupport(newProposedLinks))
            return newProposedLinks;
        else
            return proposedLinks;
    }

    private List<ProposedLinks> findNearest(List<ProposedLinks> proposedLinks)
    {
        // sorts links into shortest first and purges any conflicting longer links with conflicting breakends
        if(proposedLinks.isEmpty())
        {
            for (SvChainState svConn : mSvConnectionsMap.values())
            {
                SvVarData var = svConn.SV;

                // check whether this SV has any possible links with SVs of the same (remaining) rep count
                for (int be = SE_START; be <= SE_END; ++be)
                {
                    if (var.isSglBreakend() && be == SE_END)
                        continue;

                    boolean isStart = isStart(be);
                    double breakendPloidy = svConn.unlinked(be);

                    if (breakendPloidy == 0)
                        continue;

                    final SvBreakend breakend = var.getBreakend(isStart);

                    breakendPloidy = mLinkAllocator.getUnlinkedBreakendCount(breakend, true);

                    final List<SvLinkedPair> svLinks = mSvBreakendPossibleLinks.get(breakend);

                    if (svLinks == null)
                        continue;

                    for (final SvLinkedPair pair : svLinks)
                    {
                        if (mLinkAllocator.hasSkippedPairs(pair))
                            continue;

                        SvBreakend otherBreakend = pair.getOtherBreakend(breakend);

                        double otherBreakendPloidy = mLinkAllocator.getUnlinkedBreakendCount(otherBreakend, true);

                        if (otherBreakendPloidy == 0)
                            continue;

                        ProposedLinks proposedLink = new ProposedLinks(pair, NEAREST);
                        proposedLink.addBreakendPloidies(breakend, breakendPloidy, otherBreakend, otherBreakendPloidy);
                        proposedLinks.add(proposedLink);
                    }
                }
            }
        }

        List<ProposedLinks> shortestLinks = Lists.newArrayList();

        for(final ProposedLinks proposedLink : proposedLinks)
        {
            int index = 0;
            boolean addNew = true;
            while(index < shortestLinks.size())
            {
                final ProposedLinks otherLink = shortestLinks.get(index);

                // for proposed links with any breakend clash, just keep the shortest
                if(hasLinkClash(proposedLink.Links, otherLink.Links))
                {
                    if(proposedLink.shortestLinkDistance() < otherLink.shortestLinkDistance())
                    {
                        shortestLinks.remove(index);
                    }
                    else
                    {
                        addNew = false;
                        break;
                    }
                }
                else
                {
                    ++index;
                }
            }

            if(!addNew)
                continue;

            // insert by shortest distance first
            index = 0;
            while(index < shortestLinks.size())
            {
                if(proposedLink.shortestLinkDistance() < shortestLinks.get(index).shortestLinkDistance())
                    break;

                ++index;
            }

            LOGGER.trace("adding shortest proposed link: {} index({})", proposedLink.toString(), index);

            shortestLinks.add(index, proposedLink);
        }

        if(shortestLinks.size() > 1)
        {
            LOGGER.trace("found {} shortest non-clashing proposed links", shortestLinks.size());
        }

        return shortestLinks;
    }

    private boolean anyLinksHavePloidySupport(final List<ProposedLinks> proposedLinks)
    {
        return proposedLinks.stream().anyMatch(x -> x.hasRule(AP_SUPPORT));
    }

    private boolean checkClusterPloidySupport(final List<ProposedLinks> proposedLinks)
    {
        boolean anyLinksHasClusterPloidySupport = false;
        for(ProposedLinks proposedLink : proposedLinks)
        {
            if(proposedLink.hasRule(AP_SUPPORT))
            {
                anyLinksHasClusterPloidySupport = true;
                continue;
            }

            if(proposedLink.Links.stream().anyMatch(x -> !mPloidyLimits.linkHasPloidySupport(x, proposedLink.ploidy())))
                continue;

            proposedLink.addRule(AP_SUPPORT);
            anyLinksHasClusterPloidySupport = true;
        }

        return anyLinksHasClusterPloidySupport;
    }

    private static List<ProposedLinks> restrictProposedLinks(
            List<ProposedLinks> proposedLinks, List<ProposedLinks> newProposedLinks, ChainingRule newRule)
    {
        // now check for a match between any previously identified proposed links and this set
        if(proposedLinks.isEmpty())
            return newProposedLinks;

        if(newProposedLinks.isEmpty())
            return proposedLinks;

        for(ProposedLinks proposedLink : proposedLinks)
        {
            for (ProposedLinks newProposedLink : newProposedLinks)
            {
                if(anyLinkMatch(proposedLink.Links, newProposedLink.Links))
                {
                    proposedLink.addRule(newRule);
                    break;
                }
            }
        }

        return proposedLinks;
    }

    private static boolean anyLinkMatch(final List<SvLinkedPair> links1, final List<SvLinkedPair> links2)
    {
        for(final SvLinkedPair pair : links1)
        {
            if(links2.contains(pair))
                return true;
        }

        return false;
    }

    private void cullByPriority(List<ProposedLinks> proposedLinks)
    {
        if(proposedLinks.size() <= 1)
            return;

        // find the highest priority and remove any entries less than this
        int maxPriority = proposedLinks.stream().mapToInt(x -> x.priority()).max().getAsInt();

        int index = 0;
        while(index < proposedLinks.size())
        {
            if(proposedLinks.get(index).priority() < maxPriority)
                proposedLinks.remove(index);
            else
                ++index;
        }
    }

}
