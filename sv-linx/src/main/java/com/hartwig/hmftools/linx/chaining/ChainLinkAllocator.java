package com.hartwig.hmftools.linx.chaining;

import static java.lang.Math.max;
import static java.lang.Math.min;
import static java.lang.Math.sqrt;

import static com.hartwig.hmftools.common.variant.structural.StructuralVariantType.DUP;
import static com.hartwig.hmftools.linx.analysis.SvUtilities.copyNumbersEqual;
import static com.hartwig.hmftools.linx.analysis.SvUtilities.formatPloidy;
import static com.hartwig.hmftools.linx.chaining.ChainFinder.MIN_CHAINING_PLOIDY_LEVEL;
import static com.hartwig.hmftools.linx.chaining.ChainPloidyLimits.calcPloidyUncertainty;
import static com.hartwig.hmftools.linx.chaining.ChainPloidyLimits.ploidyOverlap;
import static com.hartwig.hmftools.linx.chaining.ChainingRule.ASSEMBLY;
import static com.hartwig.hmftools.linx.chaining.ChainingRule.FOLDBACK_SPLIT;
import static com.hartwig.hmftools.linx.chaining.ProposedLinks.CONN_TYPE_FOLDBACK;
import static com.hartwig.hmftools.linx.chaining.SvChain.reconcileChains;
import static com.hartwig.hmftools.linx.types.SvVarData.SE_END;
import static com.hartwig.hmftools.linx.types.SvVarData.SE_START;
import static com.hartwig.hmftools.linx.types.SvVarData.isStart;

import java.util.List;
import java.util.Map;

import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import com.hartwig.hmftools.linx.cn.PloidyCalcData;
import com.hartwig.hmftools.linx.types.SvBreakend;
import com.hartwig.hmftools.linx.types.SvLinkedPair;
import com.hartwig.hmftools.linx.types.SvVarData;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;

public class ChainLinkAllocator
{
    private static final Logger LOGGER = LogManager.getLogger(ChainLinkAllocator.class);

    private int mClusterId;

    private final List<SvLinkedPair> mSkippedPairs;
    private int mLinkIndex; // incrementing value for each link added to any chain
    private boolean mIsValid;
    private boolean mPairSkipped; // keep track of any excluded pair or SV without exiting the chaining routine
    private final List<SvLinkedPair> mUniquePairs; // cache of unique pairs added through c
    private int mNextChainId;
    private final Map<SvVarData, SvChainState> mSvConnectionsMap;
    private final List<SvChainState> mSvCompletedConnections;

    // references
    private final List<SvVarData> mFoldbacks;
    private final Map<SvVarData,List<SvLinkedPair>> mComplexDupCandidates;
    private final List<SvChain> mChains;
    private final Map<SvBreakend, List<SvLinkedPair>> mSvBreakendPossibleLinks;
    private final List<SvVarData> mDoubleMinuteSVs;

    public ChainLinkAllocator(
            final Map<SvBreakend, List<SvLinkedPair>> svBreakendPossibleLinks,
            final List<SvChain> chains,
            final List<SvVarData> foldbacks,
            final Map<SvVarData,List<SvLinkedPair>> complexDupCandidates,
            final List<SvVarData> doubleMinuteSVs)
    {
        mSvBreakendPossibleLinks = svBreakendPossibleLinks;
        mFoldbacks = foldbacks;
        mComplexDupCandidates = complexDupCandidates;
        mChains = chains;
        mDoubleMinuteSVs = doubleMinuteSVs;

        mSvConnectionsMap = Maps.newHashMap();
        mSvCompletedConnections = Lists.newArrayList();
        mUniquePairs = Lists.newArrayList();
        mSkippedPairs = Lists.newArrayList();
        mIsValid = true;
        mNextChainId = 0;
    }

    public final Map<SvVarData, SvChainState> getSvConnectionsMap() { return mSvConnectionsMap; }
    public final List<SvChainState> getSvCompletedConnections() { return mSvCompletedConnections; }

    public boolean hasSkippedPairs(final SvLinkedPair pair) { return mSkippedPairs.contains(pair); }
    public final List<SvLinkedPair> getUniquePairs() { return mUniquePairs; }

    public int getNextChainId() { return mNextChainId; }
    public int getLinkIndex() { return mLinkIndex; }

    public boolean isValid() { return mIsValid; }

    public boolean pairSkipped() { return mPairSkipped; }
    public void clearPairSkipped() { mPairSkipped = false; }

    public void initialise(int clusterId)
    {
        mClusterId = clusterId;

        mIsValid = true;
        mLinkIndex = 0;
        mPairSkipped = false;
        mNextChainId = 0;

        mUniquePairs.clear();
        mSkippedPairs.clear();
        mSvConnectionsMap.clear();
        mSvCompletedConnections.clear();
    }

    public void populateSvPloidyMap(final List<SvVarData> svList, boolean clusterHasReplication)
    {
        // make a cache of all unchained breakends in those of replicated SVs
        for(final SvVarData var : svList)
        {
            if(var.ploidy() <= MIN_CHAINING_PLOIDY_LEVEL)
            {
                // for now skip these
                LOGGER.debug("cluster({}) skipping SV({}) with low ploidy({} min={} max={})",
                        mClusterId, var.id(), formatPloidy(var.ploidy()), formatPloidy(var.ploidyMin()), formatPloidy(var.ploidyMax()));
                continue;
            }

            mSvConnectionsMap.put(var, new SvChainState(var, !clusterHasReplication));
        }
    }

    public void processProposedLinks(List<ProposedLinks> proposedLinksList)
    {
        boolean linkAdded = false;

        while (!proposedLinksList.isEmpty())
        {
            ProposedLinks proposedLinks = proposedLinksList.get(0);

            // in case an earlier link has invalidated the chain
            if (proposedLinks.targetChain() != null && !mChains.contains(proposedLinks.targetChain()))
                break;

            proposedLinksList.remove(0);

            if (!proposedLinks.isValid())
            {
                LOGGER.error("cluster({}) skipping invalid proposed links: {}", mClusterId, proposedLinks.toString());
                continue;
            }

            linkAdded |= addLinks(proposedLinks);

            if (!mIsValid)
                return;

            if (proposedLinks.multiConnection()) // stop after the first complex link is made
                break;
        }

        if (linkAdded)
        {
            mSkippedPairs.clear(); // any skipped links can now be re-evaluated
        }
    }

    // protected static int SPEC_LINK_INDEX = -1;
    protected static int SPEC_LINK_INDEX = 32;

    public boolean addLinks(final ProposedLinks proposedLinks)
    {
        // if a chain is specified, add the links to it
        // otherwise look for a chain which can link in these new pairs
        // and if none can be found, create a new chain with them

        // if no chain has a ploidy matching that of the new link and the new link is lower, then split the chain
        // if the chain has a lower ploidy, then only assign the ploidy of the chain
        // if the chain has a matching ploidy then recalculate it with the new SV's ploidy and uncertainty

        SvLinkedPair newPair = proposedLinks.Links.get(0);

        if (mLinkIndex == SPEC_LINK_INDEX)
        {
            LOGGER.debug("specific index({})", mLinkIndex);
        }

        final String topRule = proposedLinks.topRule().toString();
        proposedLinks.Links.forEach(x -> x.setLinkReason(topRule, mLinkIndex));

        boolean pairLinkedOnFirst = false;
        boolean addToStart = false;
        boolean linkClosesChain = false;
        boolean matchesChainPloidy = false;
        double newSvPloidy = 0;

        SvChain targetChain = null;

        if (proposedLinks.targetChain() != null)
        {
            targetChain = proposedLinks.targetChain();
            matchesChainPloidy = proposedLinks.linkPloidyMatch();
            newSvPloidy = proposedLinks.ploidy();
        }
        else if (proposedLinks.multiConnection())
        {
            // if no chain has been specified then don't search for one - this is managed by the specific rule
        }
        else
        {
            for (SvChain chain : mChains)
            {
                boolean canAddToStart = chain.canAddLinkedPair(newPair, true, true);
                boolean canAddToEnd = chain.canAddLinkedPair(newPair, false, true);

                if (!canAddToStart && !canAddToEnd)
                    continue;

                boolean couldCloseChain = (canAddToStart && canAddToEnd) ? chain.linkWouldCloseChain(newPair) : false;

                if (couldCloseChain)
                {
                    if (isDoubleMinuteDup() && mSvConnectionsMap.size() == 1 && mSvConnectionsMap.get(mDoubleMinuteSVs.get(0)) != null)
                    {
                        // allow the chain to be closed if this is the last pair other than excess DM DUP replicated SVs

                    }
                    else
                    {
                        LOGGER.trace("skipping linked pair({}) would close existing chain({})",
                                newPair.toString(), chain.id());

                        if (!mSkippedPairs.contains(newPair))
                        {
                            mPairSkipped = true;
                            mSkippedPairs.add(newPair);
                        }

                        linkClosesChain = true;
                        continue;
                    }
                }

                boolean linkedOnFirst = newPair.first() == chain.getFirstSV() || newPair.first() == chain.getLastSV();

                // final SvBreakend chainBreakend = linkedOnFirst ? newPair.firstBreakend() : newPair.secondBreakend();
                final SvBreakend newBreakend = linkedOnFirst ? newPair.secondBreakend() : newPair.firstBreakend();

                boolean ploidyMatched = copyNumbersEqual(proposedLinks.ploidy(), chain.ploidy())
                        || ploidyOverlap(proposedLinks.ploidy(), chain.ploidyUncertainty(), chain.ploidy(), chain.ploidyUncertainty());

                // check whether a match was expected
                if (!ploidyMatched)
                {
                    if (proposedLinks.linkPloidyMatch())
                        continue;

                    if (targetChain != null && targetChain.ploidy() > chain.ploidy())
                        continue; // stick with the larger ploidy chain
                }

                targetChain = chain;
                addToStart = canAddToStart;
                pairLinkedOnFirst = linkedOnFirst;

                // record the ploidy of the SV which would be added to this chain to determine whether the chain will need splitting
                newSvPloidy = proposedLinks.breakendPloidy(newBreakend);

                if (ploidyMatched)
                {
                    matchesChainPloidy = true;
                    break;
                }
            }
        }

        boolean isNewChain = (targetChain == null);
        boolean reconcileChains = !isNewChain;

        if (!isNewChain)
        {
            // scenarios:
            // - ploidy matches - add the new link and recalculate the chain ploidy
            // - foldback or complex dup with 2-1 ploidy match - replicate the chain accordingly and halve the chain ploidy
            // - foldback or complex dup with chain greater than 2x the foldback or complex dup
            //      - split off the excess and then replicate and halve the remainder
            // - foldback where the foldback itself is a chain, connecting to a single other breakend which may also be chained
            //      - here it is not the foldback chain which needs replicating or splitting but the other one
            // - normal link with chain ploidy higher - split off the chain
            // - normal link with chain ploidy lower - only allocate the chain ploidy for the new link
            boolean requiresChainSplit = false;

            if (!matchesChainPloidy && targetChain.ploidy() > newSvPloidy)
            {
                if (!proposedLinks.multiConnection())
                    requiresChainSplit = true;
                else
                    requiresChainSplit = (targetChain.ploidy() > newSvPloidy * 2);

                matchesChainPloidy = true;
            }

            if (requiresChainSplit)
            {
                SvChain newChain = new SvChain(mNextChainId++);
                mChains.add(newChain);

                // copy the existing links into a new chain and set to the ploidy difference
                newChain.copyFrom(targetChain);

                if (!proposedLinks.multiConnection())
                {
                    newChain.setPloidyData(targetChain.ploidy() - newSvPloidy, targetChain.ploidyUncertainty());
                    targetChain.setPloidyData(newSvPloidy, targetChain.ploidyUncertainty());
                }
                else if (targetChain.ploidy() > newSvPloidy * 2)
                {
                    // chain will have its ploidy halved a  nyway so just split off the excess
                    newChain.setPloidyData(targetChain.ploidy() - newSvPloidy * 2, targetChain.ploidyUncertainty());
                    targetChain.setPloidyData(newSvPloidy * 2, targetChain.ploidyUncertainty());
                }

                LOGGER.debug("new chain({}) ploidy({}) from chain({}) ploidy({}) from new SV ploidy({})",
                        newChain.id(), formatPloidy(newChain.ploidy()),
                        targetChain.id(), formatPloidy(targetChain.ploidy()), formatPloidy(newSvPloidy));
            }

            if (proposedLinks.multiConnection())
            {
                LOGGER.debug("duplicating chain({}) for multi-connect {}", targetChain.id(), proposedLinks.getSplittingRule());

                if (proposedLinks.getSplittingRule() == FOLDBACK_SPLIT)
                    targetChain.foldbackChainOnLink(proposedLinks.Links.get(0), proposedLinks.Links.get(1));
                else
                    targetChain.duplicateChainOnLink(proposedLinks.Links.get(0), proposedLinks.Links.get(1));

                double newPloidy = targetChain.ploidy() * 0.5;
                double newUncertainty = targetChain.ploidyUncertainty() * sqrt(2);
                targetChain.setPloidyData(newPloidy, newUncertainty);
            }
            else
            {
                final SvBreakend newSvBreakend = pairLinkedOnFirst ? newPair.secondBreakend() : newPair.firstBreakend();

                PloidyCalcData ploidyData;

                if (matchesChainPloidy || targetChain.ploidy() > newSvPloidy)
                {
                    if (!proposedLinks.linkPloidyMatch())
                    {
                        ploidyData = calcPloidyUncertainty(
                                new PloidyCalcData(proposedLinks.ploidy(), newSvBreakend.ploidyUncertainty()),
                                new PloidyCalcData(targetChain.ploidy(), targetChain.ploidyUncertainty()));
                    }
                    else
                    {
                        ploidyData = calcPloidyUncertainty(
                                new PloidyCalcData(proposedLinks.breakendPloidy(newSvBreakend), newSvBreakend.ploidyUncertainty()),
                                new PloidyCalcData(targetChain.ploidy(), targetChain.ploidyUncertainty()));
                    }

                    targetChain.setPloidyData(ploidyData.PloidyEstimate, ploidyData.PloidyUncertainty);
                }
                else
                {
                    // ploidy of the link is higher so keep the chain ploidy unch and reduce what can be allocated from this link
                    proposedLinks.setLowerPloidy(targetChain.ploidy());
                }

                targetChain.addLink(proposedLinks.Links.get(0), addToStart);
            }
        }
        else
        {
            if (linkClosesChain)
                return false; // skip this link for now

            // where more than one links is being added, they may not be able to be added to the same chain
            // eg a chained foldback replicating another breakend - the chain reconciliation step will join them back up
            SvChain newChain = null;
            for (final SvLinkedPair pair : proposedLinks.Links)
            {
                if (newChain != null)
                {
                    if (newChain.canAddLinkedPairToStart(pair))
                    {
                        newChain.addLink(pair, true);
                    }
                    else if (newChain.canAddLinkedPairToEnd(pair))
                    {
                        newChain.addLink(pair, false);
                    }
                    else
                    {
                        newChain = null;
                        reconcileChains = true;
                    }
                }

                if (newChain == null)
                {
                    newChain = new SvChain(mNextChainId++);
                    mChains.add(newChain);
                    targetChain = newChain;

                    newChain.addLink(pair, true);

                    PloidyCalcData ploidyData;

                    if (!proposedLinks.linkPloidyMatch() || proposedLinks.multiConnection())
                    {
                        ploidyData = calcPloidyUncertainty(
                                new PloidyCalcData(proposedLinks.ploidy(), newPair.first().ploidyUncertainty()),
                                new PloidyCalcData(proposedLinks.ploidy(), newPair.second().ploidyUncertainty()));
                    }
                    else
                    {
                        ploidyData = calcPloidyUncertainty(
                                new PloidyCalcData(proposedLinks.breakendPloidy(newPair.firstBreakend()), newPair.first()
                                        .ploidyUncertainty()),
                                new PloidyCalcData(proposedLinks.breakendPloidy(newPair.secondBreakend()), newPair.second()
                                        .ploidyUncertainty()));
                    }

                    newChain.setPloidyData(ploidyData.PloidyEstimate, ploidyData.PloidyUncertainty);
                }
            }
        }

        for (SvLinkedPair pair : proposedLinks.Links)
        {
            LOGGER.debug("index({}) method({}) adding linked pair({} ploidy={}) to {} chain({}) ploidy({})",
                    mLinkIndex, topRule, pair.toString(), formatPloidy(proposedLinks.ploidy()), isNewChain ? "new" : "existing",
                    targetChain.id(), String.format("%.1f unc=%.1f", targetChain.ploidy(), targetChain.ploidyUncertainty()));
        }

        registerNewLink(proposedLinks);
        ++mLinkIndex;

        if (reconcileChains)
        {
            // now see if any partial chains can be linked
            reconcileChains(mChains);
        }

        return true;
    }

    private void registerNewLink(final ProposedLinks proposedLink)
    {
        List<SvBreakend> exhaustedBreakends = Lists.newArrayList();
        boolean canUseMaxPloidy = proposedLink.topRule() == ASSEMBLY;

        for (final SvLinkedPair newPair : proposedLink.Links)
        {
            for (int se = SE_START; se <= SE_END; ++se)
            {
                final SvBreakend breakend = newPair.getBreakend(isStart(se));

                if (exhaustedBreakends.contains(breakend))
                    continue;

                final SvBreakend otherPairBreakend = newPair.getOtherBreakend(breakend);
                final SvVarData var = breakend.getSV();

                SvChainState svConn = mSvConnectionsMap.get(var);

                if (otherPairBreakend == null || breakend == null)
                {
                    LOGGER.error("cluster({}) invalid breakend in proposed link: {}", mClusterId, proposedLink.toString());
                    mIsValid = false;
                    return;
                }

                if (svConn == null || svConn.breakendExhaustedVsMax(breakend.usesStart()))
                {
                    LOGGER.error("breakend({}) breakend already exhausted: {} with proposedLink({})",
                            breakend.toString(), svConn != null ? svConn.toString() : "null", proposedLink.toString());
                    mIsValid = false;
                    return;
                }

                svConn.addConnection(otherPairBreakend, breakend.usesStart());

                boolean breakendExhausted = proposedLink.breakendPloidyMatched(breakend);

                if (breakendExhausted)
                {
                    // this proposed link fully allocates the breakend
                    svConn.add(breakend.usesStart(), max(svConn.unlinked(breakend.usesStart()), proposedLink.ploidy()));
                }
                else
                {
                    svConn.add(breakend.usesStart(), proposedLink.ploidy());

                    breakendExhausted = canUseMaxPloidy ? svConn.breakendExhaustedVsMax(breakend.usesStart())
                            : svConn.breakendExhausted(breakend.usesStart());
                }

                if (breakendExhausted)
                    exhaustedBreakends.add(breakend);

                final SvBreakend otherSvBreakend = var.getBreakend(!breakend.usesStart());

                if (otherSvBreakend != null)
                    removeOppositeLinks(otherSvBreakend, otherPairBreakend);
            }

            // track unique pairs to avoid conflicts (eg end-to-end and start-to-start)
            if (!matchesExistingPair(newPair))
            {
                mUniquePairs.add(newPair);
            }
        }

        // clean up breakends and SVs which have been fully allocated
        for (final SvBreakend breakend : exhaustedBreakends)
        {
            final SvVarData var = breakend.getSV();

            SvChainState svConn = mSvConnectionsMap.get(var);

            if (svConn != null)
            {
                boolean otherBreakendExhausted = canUseMaxPloidy ? svConn.breakendExhaustedVsMax(!breakend.usesStart())
                        : svConn.breakendExhausted(!breakend.usesStart());

                if (otherBreakendExhausted)
                {
                    checkSvComplete(svConn);

                    if (var.isFoldback())
                    {
                        // remove if no other instances of this SV remain
                        mFoldbacks.remove(var);
                    }
                    else if (mComplexDupCandidates.get(var) != null)
                    {
                        mComplexDupCandidates.remove(var);
                    }
                }
            }

            List<SvLinkedPair> possibleLinks = mSvBreakendPossibleLinks.get(breakend);

            if (possibleLinks != null)
            {
                // since this breakend has been exhausted, remove any links which depend on it
                removePossibleLinks(possibleLinks, breakend);
            }
        }
    }

    private void removePossibleLinks(List<SvLinkedPair> possibleLinks, SvBreakend exhaustedBreakend)
    {
        if (possibleLinks == null || possibleLinks.isEmpty())
            return;

        int index = 0;
        while (index < possibleLinks.size())
        {
            SvLinkedPair possibleLink = possibleLinks.get(index);

            if (possibleLink.hasBreakend(exhaustedBreakend))
            {
                // remove this from consideration
                possibleLinks.remove(index);

                SvBreakend otherBreakend = possibleLink.getBreakend(true) == exhaustedBreakend ?
                        possibleLink.getBreakend(false) : possibleLink.getBreakend(true);

                // and remove the pair which was cached in the other breakend's possibles list
                List<SvLinkedPair> otherPossibles = mSvBreakendPossibleLinks.get(otherBreakend);

                if (otherPossibles != null)
                {
                    for (SvLinkedPair otherPair : otherPossibles)
                    {
                        if (otherPair == possibleLink)
                        {
                            otherPossibles.remove(otherPair);

                            if (otherPossibles.isEmpty())
                            {
                                // LOGGER.debug("breakend({}) has no more possible links", otherBreakend);
                                mSvBreakendPossibleLinks.remove(otherBreakend);
                            }

                            break;
                        }
                    }
                }
            }
            else
            {
                ++index;
            }
        }

        if (possibleLinks.isEmpty())
        {
            //LOGGER.debug("breakend({}) has no more possible links", origBreakend);
            mSvBreakendPossibleLinks.remove(exhaustedBreakend);
        }
    }

    private void removeOppositeLinks(final SvBreakend otherSvBreakend, final SvBreakend otherPairBreakend)
    {
        // check for an opposite pairing between these 2 SVs - need to look into other breakends' lists

        // such a link can happen for a complex dup around a single SV, so skip if any of these exist
        if (!mComplexDupCandidates.isEmpty())
            return;

        List<SvLinkedPair> otherBeLinks = mSvBreakendPossibleLinks.get(otherSvBreakend);

        if (otherBeLinks == null)
            return;

        if (otherBeLinks.isEmpty())
        {
            mSvBreakendPossibleLinks.remove(otherSvBreakend);
            return;
        }

        final SvBreakend otherOrigBreakendAlt = otherPairBreakend.getOtherBreakend();

        if (otherOrigBreakendAlt == null)
            return;

        for (SvLinkedPair pair : otherBeLinks)
        {
            if (pair.hasBreakend(otherSvBreakend) && pair.hasBreakend(otherOrigBreakendAlt))
            {
                otherBeLinks.remove(pair);

                if (otherBeLinks.isEmpty())
                    mSvBreakendPossibleLinks.remove(otherSvBreakend);

                return;
            }
        }
    }

    private void checkSvComplete(final SvChainState svConn)
    {
        if (svConn.breakendExhausted(true) && (svConn.SV.isNullBreakend() || svConn.breakendExhausted(false)))
        {
            LOGGER.trace("SV({}) both breakends exhausted", svConn.toString());
            mSvConnectionsMap.remove(svConn.SV);
            mSvCompletedConnections.add(svConn);
        }
    }

    protected double getUnlinkedBreakendCount(final SvBreakend breakend)
    {
        SvChainState svConn = mSvConnectionsMap.get(breakend.getSV());
        if (svConn == null)
            return 0;

        return !svConn.breakendExhausted(breakend.usesStart()) ? svConn.unlinked(breakend.usesStart()) : 0;
    }

    protected double getMaxUnlinkedBreakendCount(final SvBreakend breakend)
    {
        SvChainState svConn = mSvConnectionsMap.get(breakend.getSV());
        if (svConn == null)
            return 0;

        if (!svConn.breakendExhausted(breakend.usesStart()))
            return svConn.unlinked(breakend.usesStart());
        else if (!svConn.breakendExhaustedVsMax(breakend.usesStart()))
            return svConn.maxUnlinked(breakend.usesStart());
        else
            return 0;
    }

    protected double getUnlinkedCount(final SvVarData var)
    {
        SvChainState svConn = mSvConnectionsMap.get(var);
        if (svConn == null)
            return 0;

        if (svConn.breakendExhausted(true) || svConn.breakendExhausted(false))
            return 0;

        return min(svConn.unlinked(SE_START), svConn.unlinked(SE_END));
    }

    public boolean matchesExistingPair(final SvLinkedPair pair)
    {
        for(SvLinkedPair existingPair : mUniquePairs)
        {
            if(pair.matches(existingPair))
                return true;
        }

        return false;
    }

    protected boolean isDoubleMinuteDup()
    {
        return mDoubleMinuteSVs.size() == 1 && mDoubleMinuteSVs.get(0).type() == DUP;
    }

    protected boolean isDoubleMinuteDup(final SvVarData var)
    {
        return isDoubleMinuteDup() && mDoubleMinuteSVs.get(0) == var;
    }

    private void removeSkippedPairs(List<SvLinkedPair> possiblePairs)
    {
        // some pairs are temporarily unavailable for use (eg those which would close a chain)
        // to to avoid continually trying to add them, keep them out of consideration until a new links is added
        if(mSkippedPairs.isEmpty())
            return;

        int index = 0;
        while(index < possiblePairs.size())
        {
            SvLinkedPair pair = possiblePairs.get(index);

            if(mSkippedPairs.contains(pair))
                possiblePairs.remove(index);
            else
                ++index;
        }
    }

}
