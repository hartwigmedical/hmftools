package com.hartwig.hmftools.linx.chaining;

import static java.lang.Math.max;

import static com.hartwig.hmftools.common.utils.sv.StartEndIterator.SE_END;
import static com.hartwig.hmftools.common.utils.sv.StartEndIterator.SE_PAIR;
import static com.hartwig.hmftools.common.utils.sv.StartEndIterator.SE_START;
import static com.hartwig.hmftools.common.utils.sv.StartEndIterator.isStart;
import static com.hartwig.hmftools.linx.LinxConfig.LNX_LOGGER;
import static com.hartwig.hmftools.linx.analysis.SvUtilities.formatJcn;
import static com.hartwig.hmftools.linx.chaining.ChainFinder.MIN_CHAINING_JCN_LEVEL;
import static com.hartwig.hmftools.linx.chaining.ChainJcnLimits.UNCERTAINTY_SCALE_FACTOR;
import static com.hartwig.hmftools.linx.chaining.ChainJcnLimits.calcJcnUncertainty;
import static com.hartwig.hmftools.linx.chaining.ChainJcnLimits.jcnMatch;
import static com.hartwig.hmftools.linx.chaining.ChainUtils.duplicateChainOnLink;
import static com.hartwig.hmftools.linx.chaining.ChainUtils.foldbackChainOnChain;
import static com.hartwig.hmftools.linx.chaining.ChainUtils.foldbackChainOnLink;
import static com.hartwig.hmftools.linx.chaining.ChainUtils.reconcileChains;
import static com.hartwig.hmftools.linx.chaining.ChainingRule.ASSEMBLY;
import static com.hartwig.hmftools.linx.chaining.ChainingRule.FOLDBACK_SPLIT;
import static com.hartwig.hmftools.linx.chaining.LinkSkipType.CLOSING;
import static com.hartwig.hmftools.linx.chaining.LinkSkipType.JCN_MISMATCH;

import java.util.List;
import java.util.Map;
import java.util.stream.Collectors;

import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import com.hartwig.hmftools.linx.cn.JcnCalcData;
import com.hartwig.hmftools.linx.types.LinkedPair;
import com.hartwig.hmftools.linx.types.SvBreakend;
import com.hartwig.hmftools.linx.types.SvVarData;

public class ChainLinkAllocator
{
    private int mClusterId;

    private final Map<LinkedPair,LinkSkipType> mSkippedPairs;
    private int mLinkIndex; // incrementing value for each link added to any chain
    private boolean mIsValid;
    private boolean mPairSkipped; // keep track of any excluded pair or SV without exiting the chaining routine
    private boolean mChainsSplit;
    private final List<LinkedPair> mUniquePairs; // cache of unique pairs added through chaining
    private int mNextChainId;

    // chaining state for each SV
    private final SvChainConnections mSvConnections;
    private final List<ChainState> mSvCompletedConnections; // fully exhausted SVs are moved into this collection

    // references
    private final ChainJcnLimits mJcnLimits;
    private final List<SvChain> mChains;
    private final Map<SvBreakend, List<LinkedPair>> mSvBreakendPossibleLinks;
    private final List<SvVarData> mDoubleMinuteSVs;

    public ChainLinkAllocator(
            final ChainJcnLimits jcnLimits,
            final Map<SvBreakend, List<LinkedPair>> svBreakendPossibleLinks,
            final List<SvChain> chains,
            final List<SvVarData> doubleMinuteSVs)
    {
        mJcnLimits = jcnLimits;
        mSvBreakendPossibleLinks = svBreakendPossibleLinks;
        mChains = chains;
        mDoubleMinuteSVs = doubleMinuteSVs;

        mSvConnections = new SvChainConnections();
        mSvCompletedConnections = Lists.newArrayList();
        mUniquePairs = Lists.newArrayList();
        mSkippedPairs = Maps.newHashMap();
        mIsValid = true;
        mNextChainId = 0;
    }

    public final SvChainConnections getSvConnections() { return mSvConnections; }
    public final List<ChainState> getSvCompletedConnections() { return mSvCompletedConnections; }

    public final List<LinkedPair> getUniquePairs() { return mUniquePairs; }

    public int getNextChainId() { return mNextChainId; }
    public int getLinkIndex() { return mLinkIndex; }

    public boolean isValid() { return mIsValid; }

    public boolean pairSkipped() { return mPairSkipped; }

    public void clearSkippedState()
    {
        mPairSkipped = false;
        mChainsSplit = false;
    }

    public void initialise(int clusterId)
    {
        mClusterId = clusterId;

        mIsValid = true;
        mLinkIndex = 0;
        mPairSkipped = false;
        mChainsSplit = false;
        mNextChainId = 0;

        mUniquePairs.clear();
        mSkippedPairs.clear();
        mSvConnections.clear();
        mSvCompletedConnections.clear();
    }

    public static boolean belowJcnThreshold(final SvVarData var)
    {
        return var.jcn() <= MIN_CHAINING_JCN_LEVEL;
    }

    public void populateSvJcnMap(final List<SvVarData> svList, boolean clusterHasReplication)
    {
        // make a cache of all unchained breakends in those of replicated SVs

        Double uniformClusterJcn = null;

        if(!clusterHasReplication)
        {
            uniformClusterJcn = svList.stream().mapToDouble(x-> x.jcn()).average().getAsDouble();
        }

        for(final SvVarData var : svList)
        {
            if(belowJcnThreshold(var))
            {
                // for now skip these
                LNX_LOGGER.debug("cluster({}) skipping SV({}) with low JCN({} min={} max={})",
                        mClusterId, var.id(), formatJcn(var.jcn()), formatJcn(var.jcnMin()), formatJcn(var.jcnMax()));
                continue;
            }

            mSvConnections.add(var, new ChainState(var, uniformClusterJcn));
        }
    }

    protected void addAssemblyLinksToChains(final List<LinkedPair> assembledLinks, boolean hasReplication)
    {
        if(assembledLinks.isEmpty())
            return;

        if(!hasReplication)
        {
            for (LinkedPair pair : assembledLinks)
            {
                ProposedLinks proposedLink = new ProposedLinks(pair, ASSEMBLY);
                proposedLink.addBreakendPloidies(
                        pair.getBreakend(true), getUnlinkedBreakendCount(pair.getBreakend(true)),
                        pair.getBreakend(false), getUnlinkedBreakendCount(pair.getBreakend(false)));

                if(!proposedLink.isValid())
                {
                    LNX_LOGGER.debug("cluster({}) skipping assembled link({}) with low JCN", mClusterId, proposedLink);
                    continue;
                }

                addLinks(proposedLink);
            }

            return;
        }

        // replicate any assembly links where the JCN supports it, taking note of multiple connections between the same
        // breakend and other breakends eg if a SV has JCN 2 and 2 different assembly links, it can only link once, whereas
        // if it has JCN 2 and 1 link it should be made twice, and any higher combinations are unclear

        // first gather up all the breakends which have only one assembled link and record their JCN
        List<SvBreakend> singleLinkBreakends = Lists.newArrayList();
        List<LinkedPair> bothMultiPairs = Lists.newArrayList();
        Map<SvBreakend,Double> breakendPloidies = Maps.newHashMap();

        // a working cache reduced as links as allocated
        List<LinkedPair> assemblyLinks = Lists.newArrayList(assembledLinks);

        // identify assembly links where both breakends have only 1 option, and links these immediately
        // make note of those where both breakends have multiple options
        int index = 0;
        while(index < assemblyLinks.size())
        {
            LinkedPair pair = assemblyLinks.get(index);
            final SvBreakend firstBreakend = pair.firstBreakend();
            final SvBreakend secondBreakend = pair.secondBreakend();

            boolean firstHasSingleConn = firstBreakend.getSV().getMaxAssembledBreakend() <= 1;
            boolean secondHasSingleConn = secondBreakend.getSV().getMaxAssembledBreakend() <= 1;

            double firstJcn = getUnlinkedBreakendCount(firstBreakend);
            double secondJcn = getUnlinkedBreakendCount(secondBreakend);

            if(firstJcn == 0 || secondJcn == 0)
            {
                LNX_LOGGER.debug("cluster({}) skipping assembled pair({}) with low ploidy({} & {})",
                        mClusterId, pair, formatJcn(firstJcn), formatJcn(secondJcn));

                assemblyLinks.remove(index);
                continue;
            }

            if(firstHasSingleConn && secondHasSingleConn)
            {
                ProposedLinks proposedLink = new ProposedLinks(pair, ASSEMBLY);
                proposedLink.addBreakendPloidies(firstBreakend, firstJcn, secondBreakend, secondJcn);
                addLinks(proposedLink);

                assemblyLinks.remove(index);
                continue;
            }

            ++index;

            if(firstHasSingleConn)
                singleLinkBreakends.add(firstBreakend);
            else if(secondHasSingleConn)
                singleLinkBreakends.add(secondBreakend);

            if(!firstHasSingleConn && !secondHasSingleConn)
                bothMultiPairs.add(pair);

            breakendPloidies.put(firstBreakend, firstJcn);
            breakendPloidies.put(secondBreakend, secondJcn);
        }

        // now process those pairs where one breakend has only one assembled link
        index = 0;
        while(index < assemblyLinks.size())
        {
            LinkedPair pair = assemblyLinks.get(index);

            if(!bothMultiPairs.contains(pair))
            {
                final SvBreakend firstBreakend = pair.firstBreakend();
                final SvBreakend secondBreakend = pair.secondBreakend();

                boolean firstHasSingleConn = singleLinkBreakends.contains(firstBreakend);
                boolean secondHasSingleConn = singleLinkBreakends.contains(secondBreakend);

                double firstJcn = firstHasSingleConn ? getUnlinkedBreakendCount(firstBreakend)
                        : getMaxUnlinkedBreakendCount(firstBreakend);

                double secondJcn = secondHasSingleConn ? getUnlinkedBreakendCount(secondBreakend)
                        : getMaxUnlinkedBreakendCount(secondBreakend);

                if(firstJcn == 0 || secondJcn == 0)
                {
                    LNX_LOGGER.debug("cluster({}) pair({}) assembly links already exhausted: first({}) second({})",
                            mClusterId, pair.toString(), formatJcn(firstJcn), formatJcn(secondJcn));
                    assemblyLinks.remove(index);
                    continue;
                }

                // for the breakend which has other links to make, want to avoid indicating it has been matched
                ProposedLinks proposedLink = new ProposedLinks(pair, ASSEMBLY);
                proposedLink.addBreakendPloidies(firstBreakend, firstJcn, secondBreakend, secondJcn);

                if(!firstHasSingleConn && proposedLink.exhaustBreakend(firstBreakend))
                {
                    proposedLink.overrideBreakendJcnMatched(firstBreakend, false);
                }
                else if(!secondHasSingleConn && proposedLink.exhaustBreakend(secondBreakend))
                {
                    proposedLink.overrideBreakendJcnMatched(secondBreakend, false);
                }

                LNX_LOGGER.debug("assembly multi-sgl-conn pair({}) ploidy({}): first(ploidy={} links={}) second(ploidy={} links={})",
                        pair.toString(), formatJcn(proposedLink.jcn()),
                        formatJcn(proposedLink.breakendJcn(firstBreakend)), firstBreakend.getSV().getMaxAssembledBreakend(),
                        formatJcn(proposedLink.breakendJcn(secondBreakend)), secondBreakend.getSV().getMaxAssembledBreakend());

                addLinks(proposedLink);
                assemblyLinks.remove(index);
                continue;
            }

            ++index;
        }

        // finally process the multi-connect options, most of which will now only have a single option left, and so the ploidy is known
        index = 0;
        boolean linkedPair = true;
        int iterations = 0;
        while(index < bothMultiPairs.size() && !bothMultiPairs.isEmpty())
        {
            ++iterations;
            final LinkedPair pair = bothMultiPairs.get(index);

            final SvBreakend firstBreakend = pair.firstBreakend();
            final SvBreakend secondBreakend = pair.secondBreakend();

            int firstRemainingLinks = firstBreakend.getSV().getAssembledLinkedPairs(firstBreakend.usesStart()).stream()
                    .filter(x -> bothMultiPairs.contains(x)).collect(Collectors.toList()).size();

            int secondRemainingLinks = secondBreakend.getSV().getAssembledLinkedPairs(secondBreakend.usesStart()).stream()
                    .filter(x -> bothMultiPairs.contains(x)).collect(Collectors.toList()).size();

            if(firstRemainingLinks == 0 || secondRemainingLinks == 0)
            {
                LNX_LOGGER.error("cluster({}) pair({}) unexpected remaining assembly link count: first({}) second({})",
                        mClusterId, pair.toString(), firstRemainingLinks, secondRemainingLinks);
                break;
            }

            double firstJcn = getMaxUnlinkedBreakendCount(firstBreakend);
            double secondJcn = getMaxUnlinkedBreakendCount(secondBreakend);

            if(firstJcn == 0 || secondJcn == 0)
            {
                LNX_LOGGER.debug("cluster({}) pair({}) assembly links already exhausted: first({}) second({})",
                        mClusterId, pair.toString(), formatJcn(firstJcn), formatJcn(secondJcn));
                bothMultiPairs.remove(index);
                continue;
            }

            ProposedLinks proposedLink = new ProposedLinks(pair, ASSEMBLY);

            if(firstRemainingLinks == 1 || secondRemainingLinks == 1)
            {
                proposedLink.addBreakendPloidies(firstBreakend, firstJcn, secondBreakend, secondJcn);
            }
            else if(!linkedPair)
            {
                proposedLink.addBreakendPloidies(
                        firstBreakend, firstJcn/firstRemainingLinks,
                        secondBreakend, secondJcn/secondRemainingLinks);
            }
            else
            {
                ++index;

                if(index >= bothMultiPairs.size())
                    index = 0;

                if(iterations > bothMultiPairs.size() * 3)
                {
                    LNX_LOGGER.debug("cluster({}) assembly multi-connection breakends missed", mClusterId);
                    break;
                }

                continue;
            }

            if(firstRemainingLinks > 1 && proposedLink.exhaustBreakend(firstBreakend))
            {
                proposedLink.overrideBreakendJcnMatched(firstBreakend, false);
            }

            if(secondRemainingLinks > 1 && proposedLink.exhaustBreakend(secondBreakend))
            {
                proposedLink.overrideBreakendJcnMatched(secondBreakend, false);
            }

            LNX_LOGGER.debug("assembly multi-conn pair({}) ploidy({}): first(ploidy={} links={}) second(ploidy={} links={})",
                    pair.toString(), formatJcn(proposedLink.jcn()),
                    formatJcn(firstJcn), firstBreakend.getSV().getMaxAssembledBreakend(),
                    formatJcn(secondJcn), secondBreakend.getSV().getMaxAssembledBreakend());

            addLinks(proposedLink);
            linkedPair = true;
            bothMultiPairs.remove(index);
        }

        if(!mChains.isEmpty())
        {
            LNX_LOGGER.debug("created {} partial chains from {} assembly links", mChains.size(), assembledLinks.size());
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
                LNX_LOGGER.error("cluster({}) skipping invalid proposed links: {}", mClusterId, proposedLinks.toString());
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
            if(mChainsSplit)
            {
                mSkippedPairs.clear(); // any skipped links can now be re-evaluated
            }
            else
            {
                List<LinkedPair> pairsToRemove = mSkippedPairs.entrySet().stream()
                        .filter(x -> x.getValue() != JCN_MISMATCH)
                        .map(x -> x.getKey())
                        .collect(Collectors.toList());

                pairsToRemove.stream().forEach(x -> mSkippedPairs.remove(x));
            }
        }
    }

    public boolean addLinks(final ProposedLinks proposedLinks)
    {
        // if a chain is specified, add the links to it
        // otherwise look for a chain which can link in these new pairs
        // and if none can be found, create a new chain with them

        // if no chain has a ploidy matching that of the new link and the new link is lower, then split the chain
        // if the chain has a lower ploidy, then only assign the ploidy of the chain
        // if the chain has a matching ploidy then recalculate it with the new SV's ploidy and uncertainty

        LinkedPair newPair = proposedLinks.Links.get(0);

        final String topRule = proposedLinks.topRule().toString();
        proposedLinks.Links.forEach(x -> x.setLinkReason(topRule, mLinkIndex));

        boolean addLinksToNewChain = true;
        boolean reconcileChains = false;

        if (proposedLinks.targetChain() != null)
        {
            addLinksToNewChain = false;
            addComplexLinksToExistingChain(proposedLinks);
            reconcileChains = true;
        }
        else if (proposedLinks.multiConnection())
        {
            // if no chain has been specified then don't search for one - this is managed by the specific rule-finder
        }
        else
        {
            SvChain targetChain = null;

            boolean pairLinkedOnFirst = false;
            boolean addToStart = false;
            boolean matchesChainJcn = false;
            double newSvJcn = 0;

            boolean allowChainSplits = allowDmChainSplits(proposedLinks.Links.get(0));

            boolean pairJcnMatched = ChainJcnLimits.jcnMatch(newPair.firstBreakend(), newPair.secondBreakend());

            // test every chain for whether the link would close it and look for a chain which can connect with this link

            // if one of the breakends in this new link has its other breakend in another chain and is exhausted, then force it
            // to connect to that existing chain
            // if both breakends meet this condition and the chain ploidies cannot be matched, then skip this link

            SvChain[] requiredChains = new SvChain[SE_PAIR];

            List<SvChain> allChains = Lists.newArrayList();

            for(int se = SE_START; se <= SE_END; ++se)
            {
                SvBreakend pairBreakend = newPair.getBreakend(se);

                final BreakendJcn breakendJcnData = getBreakendJcnData(pairBreakend);

                if(breakendJcnData.exhaustedInChain())
                {
                    allChains.add(breakendJcnData.MaxJcnChain);
                    requiredChains[se] = breakendJcnData.MaxJcnChain;
                    continue;
                }

                for(SvChain chain : breakendJcnData.Chains)
                {
                    if(!allChains.contains(chain))
                        allChains.add(chain);
                }
            }

            if(requiredChains[SE_START] != null || requiredChains[SE_END] != null)
            {
                if(requiredChains[SE_START] != null && requiredChains[SE_END] != null)
                {
                    if (!jcnMatch(
                            requiredChains[SE_START].jcn(), requiredChains[SE_START].jcnUncertainty(),
                            requiredChains[SE_END].jcn(), requiredChains[SE_END].jcnUncertainty()) && !pairJcnMatched)
                    {
                        if(!allowChainSplits)
                        {
                            LNX_LOGGER.trace("skipping linked pair({}) with 2 required chains({} & {}) with diff ploidies",
                                    newPair.toString(), requiredChains[SE_START].id(), requiredChains[SE_END].id());
                            addSkippedPair(newPair, JCN_MISMATCH);
                            return false;
                        }
                    }
                    else
                    {
                        matchesChainJcn = true;
                    }

                    if (requiredChains[SE_START] == requiredChains[SE_END])
                    {
                        LNX_LOGGER.trace("skipping linked pair({}) would close existing chain({})",
                                newPair.toString(), requiredChains[SE_START].id());
                        addSkippedPair(newPair, CLOSING);
                        return false;
                    }
                }

                // if both breakends are in exhausted chains, take either
                targetChain = requiredChains[SE_START] != null ? requiredChains[SE_START] : requiredChains[SE_END];

                if(targetChain.getOpenBreakend(true) == newPair.firstBreakend())
                {
                    pairLinkedOnFirst = true;
                    addToStart = true;
                }
                else if(targetChain.getOpenBreakend(true) == newPair.secondBreakend())
                {
                    pairLinkedOnFirst = false;
                    addToStart = true;
                }
                else if(targetChain.getOpenBreakend(false) == newPair.firstBreakend())
                {
                    pairLinkedOnFirst = true;
                    addToStart = false;
                }
                else if(targetChain.getOpenBreakend(false) == newPair.secondBreakend())
                {
                    pairLinkedOnFirst = false;
                    addToStart = false;
                }

                double newUncertainty = pairLinkedOnFirst ? newPair.second().jcnUncertainty() : newPair.first().jcnUncertainty();

                // check new link matches the target chain
                if(!jcnMatch(targetChain.jcn(), targetChain.jcnUncertainty(), proposedLinks.jcn(), newUncertainty)
                && !pairJcnMatched)
                {
                    if(!allowChainSplits)
                    {
                        LNX_LOGGER.trace("skipping targetChain({} ploidy={}) for proposedLink({}) on ploidy mismatch",
                                targetChain.id(), formatJcn(targetChain.jcn()), proposedLinks);

                        addSkippedPair(newPair, JCN_MISMATCH);
                        return false;
                    }
                    else
                    {
                        matchesChainJcn = false;
                        newSvJcn = proposedLinks.jcn();
                    }
                }
                else
                {
                    matchesChainJcn = true;
                    newSvJcn = targetChain.jcn();
                }

                LNX_LOGGER.trace("pair({}) links {} breakend to chain({}) as only possible connection",
                        newPair, pairLinkedOnFirst ? "first" : "second", targetChain.id());

                if(matchesChainJcn)
                {
                    // mark this breakend as to-be-exhausted unless its on both ends of the chain
                    for (int se = SE_START; se <= SE_END; ++se)
                    {
                        if (requiredChains[se] != null)
                        {
                            final SvChain requiredChain = requiredChains[se];
                            final SvBreakend breakend = newPair.getBreakend(se);

                            if (!(requiredChain.getOpenBreakend(true) == breakend && requiredChain.getOpenBreakend(false) == breakend))
                            {
                                proposedLinks.overrideBreakendJcnMatched(breakend, true);
                            }
                        }
                    }
                }
            }
            else
            {
                // look for any chain which can take this link with matching ploidy
                for (SvChain chain : allChains)
                {
                    boolean[] canAddToStart = { false, false };
                    boolean linksToFirst = false;

                    boolean ploidyMatched = jcnMatch(
                            proposedLinks.jcn(), chain.jcnUncertainty(), chain.jcn(), chain.jcnUncertainty());

                    for (int se = SE_START; se <= SE_END; ++se)
                    {
                        final SvBreakend chainBreakend = chain.getOpenBreakend(isStart(se));

                        if (chainBreakend == null)
                            continue;

                        if (chainBreakend == newPair.firstBreakend())
                            linksToFirst = true;
                        else if (chainBreakend == newPair.secondBreakend())
                            linksToFirst = false;
                        else
                            continue;

                        canAddToStart[se] = true;
                   }

                    if (!canAddToStart[SE_START] && !canAddToStart[SE_END])
                        continue;

                    boolean couldCloseChain =
                            (canAddToStart[SE_START] && canAddToStart[SE_END]) ? chain.linkWouldCloseChain(newPair) : false;

                    if (couldCloseChain)
                    {
                        LNX_LOGGER.trace("skipping linked pair({}) would close existing chain({})", newPair.toString(), chain.id());
                        addSkippedPair(newPair, CLOSING);
                        return false;
                    }

                    final SvBreakend newBreakend = linksToFirst ? newPair.secondBreakend() : newPair.firstBreakend();

                    // check whether a match was expected
                    if (!ploidyMatched)
                    {
                        if (proposedLinks.linkJcnMatch())
                            continue;

                        if (targetChain != null && targetChain.jcn() > chain.jcn())
                            continue; // stick with the larger ploidy chain
                    }
                    else if (targetChain != null && matchesChainJcn)
                    {
                        // stick with existing matched chain even if there are other equivalent options
                        continue;
                    }

                    targetChain = chain;
                    addToStart = canAddToStart[SE_START];
                    pairLinkedOnFirst = linksToFirst;

                    newSvJcn = proposedLinks.breakendJcn(newBreakend);

                    if (ploidyMatched)
                    {
                        matchesChainJcn = true;
                    }
                }
            }

            // for now don't allow chains to be split, so prevent a mismatch if the chain has a higher ploidy
            if(!allowChainSplits && targetChain != null && !matchesChainJcn && targetChain.jcn() > proposedLinks.jcn())
            {
                LNX_LOGGER.debug("skipping targetChain({} ploidy={}) for proposedLink({}) on ploidy mismatch",
                        targetChain.id(), formatJcn(targetChain.jcn()), proposedLinks);

                targetChain = null;
            }

            if (targetChain != null)
            {
                addLinksToNewChain = false;
                reconcileChains = true;
                addLinksToExistingChain(proposedLinks, targetChain, addToStart, pairLinkedOnFirst, matchesChainJcn, newSvJcn);
            }
        }

        if(addLinksToNewChain)
        {
            reconcileChains = addLinksToNewChain(proposedLinks);
        }

        if (reconcileChains)
        {
            // now see if any partial chains can be linked
            reconcileChains(mChains);
        }

        // moved to after chain reconciliation so can test whether breakends in open in chains or not
        registerNewLink(proposedLinks);
        ++mLinkIndex;

        return true;
    }

    private void addComplexLinksToExistingChain(final ProposedLinks proposedLinks)
    {
        // scenarios:
        // - ploidy matches - add the new link and recalculate the chain ploidy
        // - foldback or complex dup with 2-1 ploidy match - replicate the chain accordingly and halve the chain ploidy
        // - foldback or complex dup with chain greater than 2x the foldback or complex dup
        //      - split off the excess and then replicate and halve the remainder
        // - foldback where the foldback itself is a chain, connecting to a single other breakend which may also be chained
        //      - split the non-foldback chain and add both connections
        // - complex dup
        //      - around a single SV - just add the 2 new links
        //      - around a chain - split chain if > 2x ploidy, then add the 2 new links

        SvChain targetChain = proposedLinks.targetChain();
        boolean matchesChainJcn = proposedLinks.linkJcnMatch();
        double newSvJcn = proposedLinks.jcn();

        boolean requiresNewChain = !matchesChainJcn && targetChain.jcn() > newSvJcn * 2;

        if (requiresNewChain)
        {
            SvChain newChain = new SvChain(mNextChainId++);
            mChains.add(newChain);

            // copy the existing links into a new chain and set to the ploidy difference
            newChain.copyFrom(targetChain);

            if (targetChain.jcn() > newSvJcn * 2)
            {
                // chain will have its ploidy halved anyway so just split off the excess
                newChain.setJcnData(targetChain.jcn() - newSvJcn * 2, targetChain.jcnUncertainty());
                targetChain.setJcnData(newSvJcn * 2, targetChain.jcnUncertainty());
            }

            LNX_LOGGER.debug("new chain({}) ploidy({}) from chain({}) ploidy({}) from new SV ploidy({})",
                    newChain.id(), formatJcn(newChain.jcn()),
                    targetChain.id(), formatJcn(targetChain.jcn()), formatJcn(newSvJcn));
        }

        LNX_LOGGER.debug("duplicating chain({} links={} sv={}) for multi-connect {}",
                targetChain.id(), targetChain.getLinkCount(), targetChain.getSvCount(), proposedLinks.getSplittingRule());

        if (proposedLinks.getSplittingRule() == FOLDBACK_SPLIT)
        {
            final SvChain foldbackChain = proposedLinks.foldbackChain();

            if(foldbackChain != null)
            {
                foldbackChainOnChain(targetChain, foldbackChain, proposedLinks.Links.get(0), proposedLinks.Links.get(1));
                mChains.remove(foldbackChain);
            }
            else
            {
                foldbackChainOnLink(targetChain, proposedLinks.Links.get(0), proposedLinks.Links.get(1));
            }
        }
        else
        {
            duplicateChainOnLink(targetChain, proposedLinks.Links.get(0), proposedLinks.Links.get(1));
        }

        double newJcn = targetChain.jcn() * 0.5;
        double newUncertainty = targetChain.jcnUncertainty() / UNCERTAINTY_SCALE_FACTOR;
        targetChain.setJcnData(newJcn, newUncertainty);
        mChainsSplit = true;

        for (LinkedPair pair : proposedLinks.Links)
        {
            LNX_LOGGER.debug("index({}) method({}) adding pair({} ploidy={}) to existing chain({}) ploidy({} unc={}) match({})",
                    mLinkIndex, proposedLinks.topRule(), pair.toString(), formatJcn(proposedLinks.jcn()), targetChain.id(),
                    formatJcn(targetChain.jcn()), formatJcn(targetChain.jcnUncertainty()), proposedLinks.jcnMatchType());
        }
    }

    private void addLinksToExistingChain(final ProposedLinks proposedLinks, SvChain targetChain,
            boolean addToStart, boolean pairLinkedOnFirst, boolean matchesChainJcn, double newSvJcn)
    {
        // no longer allow chain splitting if higher ploidy than the proposed link
        final LinkedPair newPair = proposedLinks.Links.get(0);
        final SvBreakend newSvBreakend = pairLinkedOnFirst ? newPair.secondBreakend() : newPair.firstBreakend();

        if(!matchesChainJcn && targetChain.jcn() > proposedLinks.jcn())
        {
            SvChain newChain = new SvChain(mNextChainId++);
            mChains.add(newChain);

            // copy the existing links into a new chain and set to the ploidy difference
            newChain.copyFrom(targetChain);

            mChainsSplit = true;

            // chain will have its ploidy halved anyway so just split off the excess
            newChain.setJcnData(targetChain.jcn() - newSvJcn, targetChain.jcnUncertainty());
            targetChain.setJcnData(newSvJcn, targetChain.jcnUncertainty());

            LNX_LOGGER.debug("new chain({}) ploidy({}) from chain({}) ploidy({}) from new SV ploidy({})",
                    newChain.id(), formatJcn(newChain.jcn()),
                    targetChain.id(), formatJcn(targetChain.jcn()), formatJcn(newSvJcn));
        }

        JcnCalcData ploidyData;

        if (matchesChainJcn || targetChain.jcn() > newSvJcn)
        {
            if (!proposedLinks.linkJcnMatch())
            {
                ploidyData = calcJcnUncertainty(
                        new JcnCalcData(proposedLinks.jcn(), newSvBreakend.jcnUncertainty()),
                        new JcnCalcData(targetChain.jcn(), targetChain.jcnUncertainty()));
            }
            else
            {
                ploidyData = calcJcnUncertainty(
                        new JcnCalcData(proposedLinks.breakendJcn(newSvBreakend), newSvBreakend.jcnUncertainty()),
                        new JcnCalcData(targetChain.jcn(), targetChain.jcnUncertainty()));
            }

            targetChain.setJcnData(ploidyData.JcnEstimate, ploidyData.JcnUncertainty);
        }
        else
        {
            // ploidy of the link is higher so keep the chain ploidy unch and reduce what can be allocated from this link
            proposedLinks.setLowerJcn(targetChain.jcn());
        }

        targetChain.addLink(proposedLinks.Links.get(0), addToStart);

        LNX_LOGGER.debug("index({}) method({}) adding pair({} ploidy={}) to existing chain({}) ploidy({} unc={}) match({})",
                mLinkIndex, proposedLinks.topRule(), newPair.toString(), formatJcn(proposedLinks.jcn()), targetChain.id(),
                formatJcn(targetChain.jcn()), formatJcn(targetChain.jcnUncertainty()), proposedLinks.jcnMatchType());
    }

    private boolean addLinksToNewChain(final ProposedLinks proposedLinks)
    {
        boolean reconcileChains = false;
        final LinkedPair newPair = proposedLinks.Links.get(0);

        // where more than one links is being added, they may not be able to be added to the same chain
        // eg a chained foldback replicating another breakend - the chain reconciliation step will join them back up
        SvChain newChain = null;
        for (final LinkedPair pair : proposedLinks.Links)
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

                newChain.addLink(pair, true);

                JcnCalcData ploidyData;

                if (!proposedLinks.linkJcnMatch() || proposedLinks.multiConnection())
                {
                    ploidyData = calcJcnUncertainty(
                            new JcnCalcData(proposedLinks.jcn(), newPair.first().jcnUncertainty()),
                            new JcnCalcData(proposedLinks.jcn(), newPair.second().jcnUncertainty()));
                }
                else
                {
                    // blend the ploidies of the 2 SVs
                    ploidyData = calcJcnUncertainty(
                            new JcnCalcData(proposedLinks.breakendJcn(newPair.firstBreakend()),
                                    newPair.first().jcnUncertainty()),
                            new JcnCalcData(proposedLinks.breakendJcn(newPair.secondBreakend()),
                                    newPair.second().jcnUncertainty()));
                }

                newChain.setJcnData(ploidyData.JcnEstimate, ploidyData.JcnUncertainty);
            }

            LNX_LOGGER.debug("index({}) method({}) adding pair({} ploidy={}) to new chain({}) ploidy({} unc={}) match({})",
                    mLinkIndex, proposedLinks.topRule(), pair.toString(), formatJcn(proposedLinks.jcn()), newChain.id(),
                    formatJcn(newChain.jcn()), formatJcn(newChain.jcnUncertainty()), proposedLinks.jcnMatchType());
        }

        return reconcileChains;
    }

    private void registerNewLink(final ProposedLinks proposedLink)
    {
        List<SvBreakend> exhaustedBreakends = Lists.newArrayList();
        boolean canUseMaxJcn = proposedLink.topRule() == ASSEMBLY;

        for (final LinkedPair newPair : proposedLink.Links)
        {
            mJcnLimits.assignLinkJcn(newPair, proposedLink.jcn());

            mSkippedPairs.remove(newPair);

            removeOppositeLinks(newPair);

            for (int se = SE_START; se <= SE_END; ++se)
            {
                final SvBreakend breakend = newPair.getBreakend(se);

                if (exhaustedBreakends.contains(breakend))
                    continue;

                final SvBreakend otherPairBreakend = newPair.getOtherBreakend(breakend);
                final SvVarData var = breakend.getSV();

                ChainState svConn = mSvConnections.get(var);

                if (otherPairBreakend == null || breakend == null)
                {
                    LNX_LOGGER.error("cluster({}) invalid breakend in proposed link: {}", mClusterId, proposedLink.toString());
                    mIsValid = false;
                    return;
                }

                boolean beIsStart = breakend.usesStart();

                if (svConn == null || svConn.breakendExhaustedVsMax(beIsStart))
                {
                    LNX_LOGGER.error("breakend({}) breakend already exhausted: {} with proposedLink({})",
                            breakend.toString(), svConn != null ? svConn.toString() : "null", proposedLink.toString());
                    mIsValid = false;
                    return;
                }

                // scenarios:
                // 1. First connection:
                //  - ploidy-matched (most common scenario), exhaust the connection
                //  - not matched - just register ploidy against breakend
                // 2. Other breakend is exhausted:
                //  - if this breakend is now fully contained within chains, it also must be exhausted

                if(proposedLink.exhaustBreakend(breakend))
                {
                    svConn.set(beIsStart, svConn.Jcn);
                }
                else if(!svConn.hasConnections())
                {
                    // first connection
                    if(proposedLink.linkJcnMatch())
                        svConn.set(beIsStart, svConn.Jcn);
                    else
                        svConn.add(beIsStart, proposedLink.jcn());
                }
                else
                {
                    updateBreakendAllocatedPloidy(svConn, breakend);
                }

                boolean breakendExhausted = canUseMaxJcn ? svConn.breakendExhaustedVsMax(beIsStart)
                        : svConn.breakendExhausted(beIsStart);

                if (breakendExhausted)
                {
                    LNX_LOGGER.trace("{} breakend exhausted: {}", beIsStart? "start" : "end", svConn.toString());
                    exhaustedBreakends.add(breakend);
                }

                svConn.addConnection(otherPairBreakend, beIsStart);
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

            ChainState svConn = mSvConnections.get(var);

            if (svConn != null)
            {
                boolean otherBreakendExhausted = canUseMaxJcn ? svConn.breakendExhaustedVsMax(!breakend.usesStart())
                        : svConn.breakendExhausted(!breakend.usesStart());

                if (otherBreakendExhausted)
                {
                    checkSvComplete(svConn);
                }
            }

            // since this breakend has been exhausted, remove any links which depend on it
            removePossibleLinks(breakend);
        }
    }

    private void updateBreakendAllocatedPloidy(ChainState svConn, final SvBreakend breakend)
    {
        final List<SvChain> chains = mChains.stream().filter(x -> x.getSvList().contains(breakend.getSV())).collect(Collectors.toList());

        double openChainedPloidy = 0;
        double containedChainPloidy = 0;

        for(final SvChain chain : chains)
        {
            double chainPloidy = chain.jcn();

            if(chain.getOpenBreakend(true) == breakend)
                openChainedPloidy += chainPloidy;

            if(chain.getOpenBreakend(false) == breakend)
                openChainedPloidy += chainPloidy;

            containedChainPloidy += chain.getLinkedPairs().stream()
                    .filter(x -> x.hasBreakend(breakend))
                    .count() * chainPloidy;
        }

        boolean otherBreakendExhausted = svConn.breakendExhausted(!breakend.usesStart());
        if(openChainedPloidy == 0 && otherBreakendExhausted)
        {
            svConn.set(breakend.usesStart(), svConn.Jcn);
        }
        else
        {
            svConn.set(breakend.usesStart(), containedChainPloidy);
        }
    }

    private void removePossibleLinks(SvBreakend breakend)
    {
        List<LinkedPair> possibleLinks = mSvBreakendPossibleLinks.get(breakend);

        if (possibleLinks == null)
            return;

        mSvBreakendPossibleLinks.remove(breakend);

        for(LinkedPair pair : possibleLinks)
        {
            final SvBreakend otherBreakend = pair.getOtherBreakend(breakend);

            List<LinkedPair> otherPossibles = mSvBreakendPossibleLinks.get(otherBreakend);

            if (otherPossibles == null)
                continue;

            otherPossibles.remove(pair);

            if (otherPossibles.isEmpty())
                mSvBreakendPossibleLinks.remove(otherBreakend);
        }
    }

    private void removeOppositeLinks(final LinkedPair pair)
    {
        // check for an opposite pairing between these 2 SVs - need to look into other breakends' lists

        for(int se = SE_START; se <= SE_END; ++se)
        {
            SvBreakend otherBreakend = pair.getBreakend(se).getOtherBreakend();

            if(otherBreakend == null)
                continue;

            List<LinkedPair> possibleLinks = mSvBreakendPossibleLinks.get(otherBreakend);

            if (possibleLinks == null)
                continue;

            if(possibleLinks.isEmpty())
            {
                mSvBreakendPossibleLinks.remove(otherBreakend);
                continue;
            }

            SvBreakend otherPairBreakend = pair.getBreakend(!isStart(se)).getOtherBreakend();

            if(otherPairBreakend == null)
                continue;

            for (LinkedPair otherPair : possibleLinks)
            {
                if (otherPair.hasBreakend(otherBreakend) && otherPair.hasBreakend(otherPairBreakend))
                {
                    possibleLinks.remove(otherPair);

                    if (possibleLinks.isEmpty())
                        mSvBreakendPossibleLinks.remove(otherBreakend);

                    break;
                }
            }
        }
    }

    private void checkSvComplete(final ChainState svConn)
    {
        if (svConn.breakendExhausted(true) && (svConn.SV.isSglBreakend() || svConn.breakendExhausted(false)))
        {
            LNX_LOGGER.trace("SV({}) both breakends exhausted", svConn.toString());
            mSvConnections.remove(svConn.SV);
            mSvCompletedConnections.add(svConn);
        }
    }

    protected double getUnlinkedBreakendCount(final SvBreakend breakend)
    {
        ChainState svConn = mSvConnections.get(breakend.getSV());
        if (svConn == null)
            return 0;

        return !svConn.breakendExhausted(breakend.usesStart()) ? svConn.unlinked(breakend.usesStart()) : 0;
    }

    protected double getUnlinkedBreakendCount(final SvBreakend breakend, boolean limitByChains)
    {
        if(limitByChains)
            return getBreakendJcnData(breakend).unlinkedJcn();
        else
            return getUnlinkedBreakendCount(breakend);
    }

    protected BreakendJcn getBreakendJcnData(final SvBreakend breakend)
    {
        // gather up data about how much unallocated ploidy is available for this breakend
        // and whether it is tied to any chains
        ChainState svConn = mSvConnections.get(breakend.getSV());
        if (svConn == null)
            return new BreakendJcn(0, true);

        double unlinkedjcn = !svConn.breakendExhausted(breakend.usesStart()) ? svConn.unlinked(breakend.usesStart()) : 0;
        boolean otherBreakendExhausted = svConn.breakendExhausted(!breakend.usesStart());

        if(unlinkedjcn == 0)
            return new BreakendJcn(0, otherBreakendExhausted);

        if(!svConn.hasConnections())
            return new BreakendJcn(unlinkedjcn,otherBreakendExhausted);

        // if the breakend is the open end of a chain, then the chain's ploidy is a limit on the max which can assigned
        List<SvChain> chains = getChainsWithOpenBreakend(breakend);

        if(chains.isEmpty())
            return new BreakendJcn(unlinkedjcn, otherBreakendExhausted);

        // if a SV is fully connected to a chain on one side, then only offer up the chain ploidy
        if(chains.size() == 1 && otherBreakendExhausted)
        {
            final SvChain chain = chains.get(0);
            return new BreakendJcn(0, chain.jcn(), chain, chains, otherBreakendExhausted);
        }

        double totalChainjcn = 0;
        SvChain maxChain = null;
        for(final SvChain chain : chains)
        {
            if(chain.getOpenBreakend(true) == breakend)
            {
                totalChainjcn += chain.jcn();

                if(maxChain == null || maxChain.jcn() < chain.jcn())
                    maxChain = chain;
            }

            if(chain.getOpenBreakend(false) == breakend)
            {
                totalChainjcn += chain.jcn();

                if(maxChain == null || maxChain.jcn() < chain.jcn())
                    maxChain = chain;
            }
        }

        double unchainedjcn = max(unlinkedjcn - totalChainjcn, 0);

        return new BreakendJcn(unchainedjcn, totalChainjcn, maxChain, chains, false);
    }

    protected List<SvChain> getChainsWithOpenBreakend(final SvBreakend breakend)
    {
        return mChains.stream()
                .filter(x -> x.getOpenBreakend(true) == breakend || x.getOpenBreakend(false) == breakend)
                .collect(Collectors.toList());
    }

    protected double getMaxUnlinkedBreakendCount(final SvBreakend breakend)
    {
        ChainState svConn = mSvConnections.get(breakend.getSV());
        if (svConn == null)
            return 0;

        if (!svConn.breakendExhausted(breakend.usesStart()))
            return svConn.unlinked(breakend.usesStart());
        else if (!svConn.breakendExhaustedVsMax(breakend.usesStart()))
            return svConn.maxUnlinked(breakend.usesStart());
        else
            return 0;
    }

    public boolean matchesExistingPair(final LinkedPair pair)
    {
        for(LinkedPair existingPair : mUniquePairs)
        {
            if(pair.matches(existingPair))
                return true;
        }

        return false;
    }

    public boolean hasSkippedPairs(final LinkedPair pair)
    {
        return mSkippedPairs.keySet().stream().anyMatch(x -> x.matches(pair));
    }

    public int getSkippedPairCount(final LinkSkipType type)
    {
        return (int)mSkippedPairs.values().stream().filter(x -> x == type).count();
    }

    public final Map<LinkedPair, LinkSkipType> getSkippedPairs()
    {
        return mSkippedPairs;
    }

    private void addSkippedPair(final LinkedPair pair, LinkSkipType type)
    {
        if(hasSkippedPairs(pair))
            return;

        mPairSkipped = true;
        mSkippedPairs.put(pair, type);
    }

    public void removeSkippedPairs(final List<ProposedLinks> proposedLinks)
    {
        if(proposedLinks.isEmpty() || mSkippedPairs.isEmpty())
            return;

        int index = 0;
        while(index < proposedLinks.size())
        {
            if(proposedLinks.get(index).Links.stream().anyMatch(x -> hasSkippedPairs(x)))
            {
                proposedLinks.remove(index);
            }
            else
            {
                ++index;
            }
        }
    }

    public boolean allowDmChainSplits(final LinkedPair pair)
    {
        if(mDoubleMinuteSVs.isEmpty())
            return false;

        if(mDoubleMinuteSVs.contains(pair.first()) && !mDoubleMinuteSVs.contains(pair.second()))
            return true;

        if(!mDoubleMinuteSVs.contains(pair.first()) && mDoubleMinuteSVs.contains(pair.second()))
            return true;

        return false;
    }

}
