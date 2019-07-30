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

    private List<SvBreakend[]> mFoldbackBreakendPairs; // a pair of breakends forming a breakend, either chained or single
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

        mRulesToApply.add(ONLY);

        if(mHasReplication)
        {
            mRulesToApply.add(FOLDBACK_SPLIT);
            mRulesToApply.add(COMP_DUP_SPLIT);
            mRulesToApply.add(FOLDBACK);
            mRulesToApply.add(PLOIDY_MATCH);
            mRulesToApply.add(ADJACENT);
            mRulesToApply.add(PLOIDY_MAX);
        }
        else
        {
            mRulesToApply.add(ADJACENT_MATCH);
        }

        mRulesToApply.add(NEAREST);

        mFoldbacksInitialised = false;
        mFoldbackBreakendPairs.clear();
    }

    public List<ProposedLinks> findProposedLinks()
    {
        List<ProposedLinks> proposedLinks = Lists.newArrayList();

        if(mHasReplication)
            updateFoldbackBreakends();

        for(int i = 0; i < mRulesToApply.size(); ++i)
        {
            final ChainingRule rule = mRulesToApply.get(i);

            switch (rule)
            {
                case ONLY:
                    proposedLinks = findSingleOptionPairs();
                    break;

                case FOLDBACK_SPLIT:
                    proposedLinks = findFoldbackChainSplits(proposedLinks);
                    break;

                case COMP_DUP_SPLIT:
                    proposedLinks = findComplexDupSplits(proposedLinks);
                    break;

                case FOLDBACK:
                    proposedLinks = findFoldbackPairs(proposedLinks);
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

            if(i > 0 && rule != NEAREST)
                cullByPriority(proposedLinks);

            if(proposedLinks.size() == 1)
                return proposedLinks;
        }

        return proposedLinks;
    }

    private List<ProposedLinks> findSingleOptionPairs()
    {
        // find all breakends with only one other link options, or only to both breakends of the same SV - ie an INV
        List<ProposedLinks> proposedLinks = Lists.newArrayList();

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
                    /*
                    else if (proposedLink.ploidy() > otherLink.ploidy())
                    {
                        addNew = false;
                        break;
                    }
                    else
                    {
                        proposedLinks.remove(otherLink);
                    }
                    */

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

        // make allowance for a single-option foldback connection, so it can populate the rest of the rules specific info
        if(proposedLinks.size() == 1 && !mFoldbackBreakendPairs.isEmpty())
        {
            for(final SvBreakend[] breakendPair : mFoldbackBreakendPairs)
            {
                SvBreakend foldbackStart = breakendPair[SE_START];
                SvBreakend foldbackEnd = breakendPair[SE_END];

                if (foldbackEnd == foldbackStart || foldbackEnd.getSV() == foldbackStart.getSV())
                {
                    final SvLinkedPair proposedLink = proposedLinks.get(0).Links.get(0);

                    if (proposedLink.hasBreakend(foldbackStart) || proposedLink.hasBreakend(foldbackEnd))
                    {
                        proposedLinks = findFoldbackChainSplits(proposedLinks);
                        break;
                    }
                }
            }
        }

        return proposedLinks;
    }

    private void updateFoldbackBreakends()
    {
        List<SvBreakend> processedBreakends = Lists.newArrayList(); // to avoid double-counting
        List<SvBreakend[]> existingChainedPairs = Lists.newArrayList();

        if(!mFoldbacksInitialised)
        {
            mFoldbacksInitialised = true;
            mFoldbackBreakendPairs.clear();

            for (SvVarData foldback : mFoldbacks)
            {
                if (foldback.isSingleBreakendFoldback())
                    continue;

                boolean isChainedFoldback = foldback.isChainedFoldback();

                SvBreakend foldbackStart = null;
                SvBreakend foldbackEnd = null;

                if (!isChainedFoldback)
                {
                    foldbackStart = foldback.getBreakend(true);
                    foldbackEnd = foldback.getBreakend(false);
                }
                else
                {
                    if (foldback.getFoldbackBreakend(true) != null)
                    {
                        foldbackStart = foldback.getBreakend(true);
                        foldbackEnd = foldback.getFoldbackBreakend(true);
                    }
                    else
                    {
                        foldbackStart = foldback.getBreakend(false);
                        foldbackEnd = foldback.getFoldbackBreakend(false);
                    }

                    if (processedBreakends.contains(foldbackStart) || processedBreakends.contains(foldbackEnd))
                        continue;
                }

                if (mLinkAllocator.getUnlinkedBreakendCount(foldbackStart) == 0
                || mLinkAllocator.getUnlinkedBreakendCount(foldbackEnd) == 0)
                    continue;

                final SvBreakend[] breakendPair = { foldbackStart, foldbackEnd };
                mFoldbackBreakendPairs.add(breakendPair);

                processedBreakends.add(foldbackStart);
                processedBreakends.add(foldbackEnd);
            }
        }
        else
        {
            int index = 0;
            while(index < mFoldbackBreakendPairs.size())
            {
                final SvBreakend[] breakendPair = mFoldbackBreakendPairs.get(index);

                SvBreakend foldbackStart = breakendPair[SE_START];
                SvBreakend foldbackEnd = breakendPair[SE_END];

                if(mLinkAllocator.getUnlinkedBreakendCount(foldbackStart) == 0 || mLinkAllocator.getUnlinkedBreakendCount(foldbackEnd) == 0)
                {
                    LOGGER.debug("foldback breakends({} & {}) removed from consideration", foldbackStart, foldbackEnd);
                    mFoldbackBreakendPairs.remove(index);
                    continue;
                }

                // will be checked against current chains next up
                if(foldbackStart == foldbackEnd || foldbackStart.getSV() != foldbackEnd.getSV())
                {
                    existingChainedPairs.add(breakendPair);
                }

                processedBreakends.add(foldbackStart);
                processedBreakends.add(foldbackEnd);

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

            if(mLinkAllocator.getUnlinkedBreakendCount(chainStart) == 0 || mLinkAllocator.getUnlinkedBreakendCount(chainEnd) == 0)
                continue;

            // are breakends consecutive and facing the same direction?
            if(!chainStart.chromosome().equals(chainEnd.chromosome()) || chainStart.orientation() != chainEnd.orientation())
                continue;

            if(chainStart != chainEnd && abs(chainStart.getClusterChrPosIndex() - chainEnd.getClusterChrPosIndex()) != 1)
                continue;

            final SvBreakend[] breakendPair = { chainStart, chainEnd };

            removeBreakendPair(existingChainedPairs, breakendPair);

            // if(processedBreakends.contains(chainStart) || processedBreakends.contains(chainEnd))
            //    continue;

            if(!containsBreakendPair(mFoldbackBreakendPairs, breakendPair))
            {
                LOGGER.debug("chain({}) adding chained foldback breakends({} & {})", chain.id(), chainStart, chainEnd);
                mFoldbackBreakendPairs.add(breakendPair);
            }
        }

        // remove any chained foldbacks no longer supported by the current set of chains
        existingChainedPairs.stream().forEach(x -> removeBreakendPair(mFoldbackBreakendPairs, x));
    }

    private boolean containsBreakendPair(final List<SvBreakend[]> pairs, final SvBreakend[] pair)
    {
        return pairs.stream().anyMatch(x -> (x[0] == pair[0] && x[1] == pair[1]) || (x[0] == pair[1] && x[1] == pair[0]));
    }

    private void removeBreakendPair(final List<SvBreakend[]> pairs, final SvBreakend[] pair)
    {
        for(int i = 0; i < pairs.size(); ++i)
        {
            SvBreakend[] otherPair = pairs.get(i);

            if ((otherPair[0] == pair[0] && otherPair[1] == pair[1]) || (otherPair[0] == pair[1] && otherPair[1] == pair[0]))
            {
                pairs.remove(i);
                return;
            }
        }
    }

    private List<ProposedLinks> findFoldbackChainSplits(final List<ProposedLinks> proposedLinks)
    {
        updateFoldbackBreakends();

        // both ends of a foldback connect to one end of another SV with ploidy >= 2x
        // of the other breakend is in a chain, make note of this and use the chain's ploidy
        if(mFoldbackBreakendPairs.isEmpty())
            return proposedLinks;

        List<ProposedLinks> newProposedLinks = Lists.newArrayList();

        for(final SvBreakend[] breakendPair : mFoldbackBreakendPairs)
        {
            SvBreakend foldbackStart = breakendPair[SE_START];
            SvBreakend foldbackEnd = breakendPair[SE_END];

            SvVarData foldback = foldbackStart.getSV();

            // a chained foldback is now either 2 different SVs or a chain with the same breakend exposed on both ends
            boolean isSingleBreakendFoldback = foldbackStart == foldbackEnd;
            boolean isChainedFoldback = foldbackStart.getSV() != foldbackEnd.getSV() || isSingleBreakendFoldback;

            double origFoldbackPloidy = foldback.ploidy();

            double foldbackPloidy = min(
                    mLinkAllocator.getUnlinkedBreakendCount(foldbackStart),
                    mLinkAllocator.getUnlinkedBreakendCount(foldbackEnd));

            double foldbackUncertainty = foldback.ploidyUncertainty();

            if(foldbackPloidy == 0)
                continue;

            // if a foldback breakend is chained, use the ploidy of the chain
            SvChain foldbackChain = null;

            if(isChainedFoldback)
            {
                for (SvChain chain : mChains)
                {
                    final SvBreakend chainStart = chain.getOpenBreakend(true);
                    final SvBreakend chainEnd = chain.getOpenBreakend(false);

                    if ((chainStart == foldbackStart && chainEnd == foldbackEnd) || (chainStart == foldbackEnd && chainEnd == foldbackStart))
                    {
                        foldbackChain = chain;
                        break;
                    }
                }

                if(foldbackChain != null)
                {
                    foldbackPloidy = foldbackChain.ploidy();
                    foldbackUncertainty = foldbackChain.ploidyUncertainty();
                }
            }

            List<SvLinkedPair> pairsOnFbStart = mSvBreakendPossibleLinks.get(foldbackStart);
            List<SvLinkedPair> pairsOnFbEnd = mSvBreakendPossibleLinks.get(foldbackEnd);

            if(pairsOnFbStart == null || pairsOnFbEnd == null)
                continue;

            // get the set of possible pairs for the start and end breakends of the foldback and look for:
            // a) the same non-foldback breakend linked to different breakends of the same foldback
            // b) at least a 2:1 ploidy ratio between the non-foldback breakend and the foldback breakend

            for(SvLinkedPair pairStart : pairsOnFbStart)
            {
                SvVarData nonFbVar = pairStart.getOtherSV(foldback);
                SvBreakend otherBreakend = pairStart.getOtherBreakend(foldbackStart);

                if(nonFbVar.ploidy() < origFoldbackPloidy)
                    continue;

                double nonFoldbackPloidy = mLinkAllocator.getUnlinkedBreakendCount(otherBreakend);

                // foldback ploidy must be half or less to match for a potential chain split and the foldback cannot
                // be split itself so its entire ploidy must be exhausted when it operates as a foldback
                if(!ploidyMatch(foldbackPloidy * 2, foldbackUncertainty, nonFoldbackPloidy, nonFbVar.ploidyUncertainty())
                && nonFoldbackPloidy <= foldbackPloidy * 2)
                {
                    continue;
                }

                // find the other pairing - would expect this to exist
                SvLinkedPair pairEnd = pairsOnFbEnd.stream().filter(x -> x.getOtherBreakend(foldbackEnd) == otherBreakend)
                        .findFirst().orElse(null);

                if(pairEnd == null)
                    continue;

                // if this breakend is in a chain and matches 2x the foldback ploidy, then target this chain
                // if the chain exeeds the ploidy of the foldback, then look for a higher unlinked breakend

                double maxNonFoldbackChainPloidy = 0;
                double totalNonFoldbackChainPloidy = 0;
                SvChain targetChain = null;
                boolean chainPloidyMatched = false;

                for (SvChain chain : mChains)
                {
                    final SvBreakend chainStart = chain.getOpenBreakend(true);
                    final SvBreakend chainEnd = chain.getOpenBreakend(false);

                    if(chainStart == otherBreakend || chainEnd == otherBreakend)
                    {
                        if(ploidyMatch(foldbackPloidy * 2, foldbackUncertainty, chain.ploidy(), chain.ploidyUncertainty()))
                        {
                            chainPloidyMatched = true;
                            targetChain = chain;
                            break;
                        }

                        // take a chain which could be split as second-best option
                        if(chain.ploidy() > foldbackPloidy * 2)
                        {
                            targetChain = chain;
                        }

                        maxNonFoldbackChainPloidy = max(maxNonFoldbackChainPloidy, chain.ploidy());

                        if(chainStart == otherBreakend)
                            totalNonFoldbackChainPloidy += chain.ploidy();

                        if(chainEnd == otherBreakend)
                            totalNonFoldbackChainPloidy += chain.ploidy();
                    }
                }

                double linkPloidy = 0;
                if(chainPloidyMatched)
                {
                    linkPloidy = targetChain.ploidy();
                }
                else
                {
                    // work out whether the unlinked & unchained breakend ploidy exceeds the largest viable chain
                    double unchainedPloidy = max(nonFoldbackPloidy - totalNonFoldbackChainPloidy, 0);

                    // a better connection can be made to the unchained breakend for this SV
                    if(unchainedPloidy > maxNonFoldbackChainPloidy)
                    {
                        // compare ploidies again
                        if(ploidyMatch(foldbackPloidy * 2, foldbackUncertainty, unchainedPloidy, nonFbVar.ploidyUncertainty())
                        || unchainedPloidy > foldbackPloidy * 2)
                        {
                            targetChain = null;
                            linkPloidy = unchainedPloidy;
                        }
                    }

                    /*
                    if(linkPloidy == 0)
                    {
                        // can the link be made between a chain and an open end of this SV
                        if(ploidyMatch(foldbackPloidy, foldbackUncertainty, maxNonFoldbackChainPloidy, nonFbVar.ploidyUncertainty())
                        && (ploidyMatch(foldbackPloidy, foldbackUncertainty, unchainedPloidy, nonFbVar.ploidyUncertainty())
                            || unchainedPloidy > foldbackPloidy))
                        {
                            linkPloidy = maxNonFoldbackChainPloidy;
                        }
                    }
                    */
                }

                if(linkPloidy == 0)
                    continue;

                // for a simple INV foldback which can be inserted into an existing chain, then limit the ploidy to that chain
                // and specify the target chain - otherwise the links will go into new chains and rely on chain merging

                ProposedLinks proposedLink = new ProposedLinks(
                            Lists.newArrayList(pairStart, pairEnd), FOLDBACK_SPLIT, targetChain, foldbackChain);

                proposedLink.addFoldbackBreakends(
                        foldbackStart, foldbackEnd, foldbackPloidy,
                        otherBreakend, nonFoldbackPloidy, linkPloidy, otherBreakend.ploidyUncertainty());

                newProposedLinks.add(proposedLink);

                LOGGER.trace("foldback breakends({} & {}) ploidy({}) matched with breakend({}) ploidy({}) linkPloidy({})",
                        foldbackStart, foldbackEnd, formatPloidy(foldbackPloidy),
                        otherBreakend, formatPloidy(nonFoldbackPloidy), formatPloidy(linkPloidy));
            }
        }

        // now check for a match between any previously identified proposed links and this set
        // note that to preserve the additional proposed-link info for complex links, the higher rule is considered 'new'
        //return restrictProposedLinks(newProposedLinks, proposedLinks, ONLY);

        // now check for a match between any previously identified proposed links and this set
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

    private List<ProposedLinks> findFoldbackPairs(final List<ProposedLinks> proposedLinks)
    {
        // look for 2 foldbacks facing each other
        if(mFoldbackBreakendPairs.isEmpty())
            return proposedLinks;

        if(!proposedLinks.isEmpty())
        {
            for(ProposedLinks proposedLink : proposedLinks)
            {
                final SvLinkedPair pair = proposedLink.Links.get(0);

                if(mFoldbackBreakendPairs.stream().anyMatch(x -> x[0] == pair.firstBreakend() || x[1] == pair.firstBreakend())
                && mFoldbackBreakendPairs.stream().anyMatch(x -> x[0] == pair.secondBreakend() || x[1] == pair.secondBreakend()))
                {
                    proposedLink.addRule(FOLDBACK);
                }
            }

            return proposedLinks;
        }

        List<ProposedLinks> newProposedLinks = Lists.newArrayList();

        for(final SvBreakend[] breakendPair : mFoldbackBreakendPairs)
        {
            for(int se = SE_START; se <= SE_END; ++se)
            {
                SvBreakend breakend = breakendPair[se];

                double foldbackPloidy = mLinkAllocator.getUnlinkedBreakendCount(breakend, true);

                if (foldbackPloidy == 0)
                    continue;

                final List<SvLinkedPair> svLinks = mSvBreakendPossibleLinks.get(breakend);

                if (svLinks == null)
                    continue;

                for (final SvLinkedPair pair : svLinks)
                {
                    SvBreakend otherBreakend = pair.getOtherBreakend(breakend);

                    if(!mFoldbackBreakendPairs.stream()
                            .filter(x -> x != breakendPair)
                            .anyMatch(x -> x[0] == otherBreakend || x[1] == otherBreakend))
                    {
                        continue;
                    }

                    double otherBreakendPloidy = mLinkAllocator.getUnlinkedBreakendCount(otherBreakend, true);

                    if (otherBreakendPloidy == 0)
                        continue;

                    if (mLinkAllocator.hasSkippedPairs(pair))
                        continue;

                    if (proposedLinks.stream().map(x -> x.Links.get(0)).anyMatch(y -> y == pair))
                        continue;

                    ProposedLinks proposedLink = new ProposedLinks(pair, FOLDBACK);
                    proposedLink.addBreakendPloidies(breakend, foldbackPloidy, otherBreakend, otherBreakendPloidy);

                    boolean addNew = true;

                    if(addNew)
                    {
                        LOGGER.trace("pair({}) of foldbacks with ploidy({} & {})",
                                pair.toString(), formatPloidy(foldbackPloidy), formatPloidy(otherBreakendPloidy));

                        newProposedLinks.add(proposedLink);
                    }
                }
            }
        }

        // now check for a match between any previously identified proposed links and this set
        return restrictProposedLinks(proposedLinks, newProposedLinks, FOLDBACK);
    }

    private List<ProposedLinks> findPloidyMatchPairs(List<ProposedLinks> proposedLinks)
    {
        // find pairs of matching ploidy breakends, taking the shortest where multiple exist
        List<ProposedLinks> newProposedLinks = Lists.newArrayList();

        if(!proposedLinks.isEmpty())
        {
            proposedLinks.stream().filter(x -> x.ploidyMatchType() == PM_MATCHED).forEach(x -> x.addRule(PLOIDY_MATCH));
            proposedLinks.stream().filter(x -> x.ploidyMatchType() == PM_OVERLAP).forEach(x -> x.addRule(PLOIDY_OVERLAP));
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

        return newProposedLinks;
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

        if(!proposedLinks.isEmpty())
        {
            for(ProposedLinks proposedLink : proposedLinks)
            {
                // check for any links which are in the adjacent set
                if(proposedLink.Links.stream().anyMatch(x -> mAdjacentPairs.contains(x)))
                {
                    proposedLink.addRule(ADJACENT);
                }
            }

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

        return restrictProposedLinks(proposedLinks, newProposedLinks, ADJACENT);
    }

    private List<ProposedLinks> findHighestPloidy(List<ProposedLinks> proposedLinks)
    {
        List<ProposedLinks> newProposedLinks = Lists.newArrayList();

        if(!proposedLinks.isEmpty())
        {
            // take the highest from amongst the proposed links
            double maxPloidy = proposedLinks.stream().mapToDouble(x -> x.ploidy()).max().getAsDouble();

            proposedLinks.stream().filter(x -> copyNumbersEqual(maxPloidy, x.ploidy())).forEach(x -> x.addRule(PLOIDY_MAX));
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

        return newProposedLinks;
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

    private void setClusterPloidySupport(ProposedLinks proposedLinks)
    {
        if(proposedLinks.hasRule(AP_SUPPORT))
            return;

        if(proposedLinks.Links.stream().anyMatch(x -> !mPloidyLimits.linkHasPloidySupport(x, proposedLinks.ploidy())))
            return;

        proposedLinks.addRule(AP_SUPPORT);
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

        proposedLinks.stream().forEach(x -> setClusterPloidySupport(x));

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
