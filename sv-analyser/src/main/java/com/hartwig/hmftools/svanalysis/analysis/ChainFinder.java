package com.hartwig.hmftools.svanalysis.analysis;

import static java.lang.Math.abs;
import static java.lang.Math.max;
import static java.lang.Math.min;
import static java.lang.Math.round;

import static com.hartwig.hmftools.common.variant.structural.StructuralVariantType.SGL;
import static com.hartwig.hmftools.svanalysis.analysis.LinkFinder.MIN_TEMPLATED_INSERTION_LENGTH;
import static com.hartwig.hmftools.svanalysis.analysis.LinkFinder.areLinkedSection;
import static com.hartwig.hmftools.svanalysis.analysis.SvUtilities.getProximity;
import static com.hartwig.hmftools.svanalysis.types.SvVarData.SVI_END;
import static com.hartwig.hmftools.svanalysis.types.SvVarData.SVI_START;
import static com.hartwig.hmftools.svanalysis.types.SvVarData.isStart;
import static com.hartwig.hmftools.svanalysis.types.SvLinkedPair.LINK_TYPE_TI;

import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.stream.Collectors;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.svanalysis.types.SvBreakend;
import com.hartwig.hmftools.svanalysis.types.SvChain;
import com.hartwig.hmftools.svanalysis.types.SvCluster;
import com.hartwig.hmftools.svanalysis.types.SvVarData;
import com.hartwig.hmftools.svanalysis.types.SvLinkedPair;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;

public class ChainFinder
{
    private static final Logger LOGGER = LogManager.getLogger(ChainFinder.class);

    private SvCluster mCluster;

    // chaining state
    private List<SvLinkedPair> mAssemblyLinkedPairs;
    private List<SvChain> mCompleteChains;
    private List<SvChain> mIncompleteChains;

    private List<SvChain> mPartialChains;
    private List<SvVarData> mFullyLinkedSVs;
    private List<SvVarData> mUnlinkedSVs;
    private List<SvBreakend> mUnlinkedBreakends;
    private List<SvLinkedPair> mChainedPairs;
    private Map<SvVarData,Integer> mSvReplicationMap;

    private Map<SvBreakend,List<SvLinkedPair>> mSvBreakendPossibleLinks;

    private boolean mUseNewMethod;
    private int mReqChainCount;
    private boolean mLogVerbose;

    public ChainFinder()
    {
        mCompleteChains = Lists.newArrayList();
        mIncompleteChains = Lists.newArrayList();
        mPartialChains = Lists.newArrayList();
        mUnlinkedSVs = Lists.newArrayList();
        mUnlinkedBreakends = Lists.newArrayList();
        mFullyLinkedSVs = Lists.newArrayList();
        mChainedPairs = Lists.newArrayList();
        mSvReplicationMap = new HashMap();
        mSvBreakendPossibleLinks = new HashMap();

        mLogVerbose = false;
        mReqChainCount = 0;
        mUseNewMethod = false;
    }

    public void initialise(SvCluster cluster)
    {
        mCluster = cluster;
        mReqChainCount = 0;

        mAssemblyLinkedPairs = Lists.newArrayList(cluster.getAssemblyLinkedPairs());

        mCompleteChains.clear();
        mIncompleteChains.clear();
        mPartialChains.clear();
        mChainedPairs.clear();
        mFullyLinkedSVs.clear();
        mUnlinkedSVs.clear();
        mUnlinkedBreakends.clear();
    }

    public void setLogVerbose(boolean toggle) { mLogVerbose = toggle; }
    public void setUseNewMethod(boolean toggle) { mUseNewMethod = toggle; }

    public boolean formClusterChains()
    {
        // take the assembly links as a given and then try out the inferred links to see if a single chain can be formed from all the breakends

        // only factor in templated insertions to form chains, even if they also exhibit DBs
        List<SvVarData> svList = Lists.newArrayList(mCluster.getSVs());

        if(svList.size() < 2)
            return false;

        mReqChainCount = svList.size();

        if (mCluster.getCount() >= 4)
        {
            LOGGER.debug("cluster({}) assemblyLinks({}) svCount({} all={}) existingChains({})",
                    mCluster.id(), mAssemblyLinkedPairs.size(), svList.size(), mCluster.getCount(), mCluster.getChains().size());
        }

        if(mUseNewMethod)
        {
            buildSvChains();

            for(final SvChain chain : mPartialChains)
            {
                mCluster.addChain(chain);
            }

        }

        List<SvChain> chains = Lists.newArrayList();
        findSvChainsIncrementally(svList, chains);

        if(chains.size() == 1 && chains.get(0).getSvCount() == mReqChainCount)
        {
            SvChain completeChain = chains.get(0);
            mCompleteChains.add(completeChain);
        }
        else
        {
            for (SvChain chain : chains)
            {
                if (chain.getSvCount() < 2)
                    continue;

                mIncompleteChains.add(chain);
            }
        }

        // first check for any complete chains, and if found, add the shortest one
        boolean fullyChained = cacheCompleteChain();

        if(!fullyChained)
        {
            // otherwise add the longest mutually exclusive chains
            cacheIncompleteChains();
        }

        return fullyChained;
    }

    private boolean cacheCompleteChain()
    {
        if(mCompleteChains.isEmpty())
            return false;

        // cache any complete chain found
        int shortestLength = -1;
        SvChain shortestFullChain = null;
        for (SvChain chain : mCompleteChains)
        {
            chain.recalcLength();

            if (chain.getSvCount() != mReqChainCount)
                continue;

            if (shortestFullChain == null || chain.getLength() < shortestLength)
            {
                shortestFullChain = chain;
                shortestLength = chain.getLength();
            }
        }

        if (shortestFullChain == null)
            return false;

        shortestFullChain.setId(mCluster.getChains().size());

        mCluster.getChains().clear();
        mCluster.addChain(shortestFullChain);
        return true;
    }

    private void cacheIncompleteChains()
    {
        if(mIncompleteChains.isEmpty())
            return;

        mCluster.getChains().clear();

        List<SvVarData> chainedSVs = Lists.newArrayList();

        // otherwise add the longest mutually exclusive chains
        while (!mIncompleteChains.isEmpty())
        {
            int maxSvLength = -1;
            SvChain maxLengthChain = null;

            for (SvChain chain : mIncompleteChains)
            {
                // cannot contain existing chained SV
                boolean hasExistingSVs = false;

                if (!chainedSVs.isEmpty())
                {
                    for (final SvVarData var : chain.getSvList())
                    {
                        if (chainedSVs.contains(var))
                        {
                            hasExistingSVs = true;
                            break;
                        }
                    }
                }

                if (hasExistingSVs)
                    continue;

                if (chain.getSvCount() > maxSvLength)
                {
                    maxSvLength = chain.getSvCount();
                    maxLengthChain = chain;
                }
            }

            if (maxLengthChain == null)
                break;

            maxLengthChain.setId(mCluster.getChains().size());

            if(maxLengthChain.getLinkCount() >= 3 || mLogVerbose)
            {
                LOGGER.debug("cluster({}) adding incomplete chain({}) length({}) with {} linked pairs",
                        mCluster.id(), maxLengthChain.id(), maxLengthChain.getLength(), maxLengthChain.getLinkCount());

                maxLengthChain.logLinks();
            }

            mCluster.addChain(maxLengthChain);
            chainedSVs.addAll(maxLengthChain.getSvList());
            mIncompleteChains.remove(maxLengthChain);
        }
    }

    private void buildSvChains()
    {
        setUnlinkedBreakends();

        // first make chains out of any assembly links
        addAssemblyLinksToChains();

        determinePossibleLinks();

        setSvReplicationCounts();

        // now find the best next candidate link giving priority to replication count, following by min link resolution
        // and lastly shortest distance

        // after each new link, attempt to link any partial chains together
        boolean isFinished = false;
        while(!isFinished)
        {
            List<SvLinkedPair> possiblePairs = null;

            if (!mSvReplicationMap.isEmpty())
            {
                // find highest replication count if there is one
                List<SvVarData> maxRepSVs = getMaxReplicationSvIds();

                if (!maxRepSVs.isEmpty())
                {
                    List<SvBreakend> breakendList = mUnlinkedBreakends.stream()
                            .filter(x -> !x.getSV().isReplicatedSv())
                            .filter(x -> maxRepSVs.contains(x.getOrigSV()))
                            .collect(Collectors.toList());

                    possiblePairs = findRestrictedPairs(breakendList, true);
                    // possiblePairs = findRestrictedPairs(maxRepSVs, true);
                }
            }

            if (possiblePairs == null)
            {
                List<SvBreakend> breakendList = mUnlinkedBreakends.stream()
                        .filter(x -> !x.getSV().isReplicatedSv())
                        .collect(Collectors.toList());

                possiblePairs = findRestrictedPairs(breakendList, false);
                // possiblePairs = findRestrictedPairs(mUnlinkedSVs, false);
            }

            if (possiblePairs.isEmpty())
            {
                isFinished = true;
                break;
            }

            // add the shortest of these - or all which are mutually exclusive??
            while (!possiblePairs.isEmpty())
            {
                SvLinkedPair shortestPair = null;
                for (SvLinkedPair pair : possiblePairs)
                {
                    if (shortestPair == null || pair.length() < shortestPair.length())
                    {
                        shortestPair = pair;
                    }
                }

                if (shortestPair == null)
                {
                    isFinished = true;
                    break;
                }

                possiblePairs.remove(shortestPair);

                // remove any other conflicting pairs
                int index = 0;
                while (index < possiblePairs.size())
                {
                    SvLinkedPair pair = possiblePairs.get(index);
                    if (pair.hasLinkClash(shortestPair) || shortestPair.oppositeMatch(pair))
                    {
                        possiblePairs.remove(index);
                    }
                    else
                    {
                        ++index;
                    }
                }

                // the actual pair being added by be drawn from unlinked breakends, not the real ones
                addPairToChain(shortestPair, false);
            }
        }

        if(mLogVerbose)
        {
            LOGGER.debug("cluster({}) chaining finished: chains({}) unlinked(SVs={} breakends={})",
                    mCluster.id(), mPartialChains.size(), mUnlinkedSVs.size(), mUnlinkedBreakends.size());

            for(final SvChain chain : mPartialChains)
            {
                LOGGER.debug("cluster({}) adding chain({}) with {} linked pairs:",
                        mCluster.id(), chain.id(), chain.getLinkCount());
                chain.logLinks();
            }
        }
    }

    private List<SvLinkedPair> findRestrictedPairs(List<SvBreakend> breakendList, boolean isRestricted)
    {
        // of these pairs, do some have less alternatives links which could be made than others
        // eg if high-rep SVs are A and B, and possible links have been found A-C, A-D and B-E,
        // then count how many links could be made between C, D and E to other SVs, and then select the least restrictive
        // say A-C is the only link C can make and B-D is the only link that D can make, then would want to make them both
        int minPairCount = 0;
        List<SvLinkedPair> minLinkPairs = Lists.newArrayList();

        for(SvBreakend breakend : breakendList)
        {
            List<SvLinkedPair> possiblePairs = mSvBreakendPossibleLinks.get(breakend);

            if (possiblePairs == null || possiblePairs.isEmpty())
                continue;

            if (minPairCount == 0 || possiblePairs.size() < minPairCount)
            {
                minLinkPairs.clear();
                minPairCount = possiblePairs.size();
            }

            if (possiblePairs.size() == minPairCount)
            {
                for (SvLinkedPair pair : possiblePairs)
                {
                    if (!minLinkPairs.contains(pair))
                        minLinkPairs.add(pair);
                }
            }

            if (isRestricted)
            {
                // also check this max-rep SV's pairings to test how they are restricted
                for (SvLinkedPair pair : possiblePairs)
                {
                    SvBreakend otherBreakend = pair.getOtherBreakend(breakend);

                    List<SvLinkedPair> otherPossiblePairs = mSvBreakendPossibleLinks.get(otherBreakend);

                    if (otherPossiblePairs == null || otherPossiblePairs.isEmpty())
                        continue; // logical assert

                    if (minPairCount == 0 || otherPossiblePairs.size() < minPairCount)
                    {
                        minLinkPairs.clear();
                        minPairCount = otherPossiblePairs.size();
                    }

                    if (otherPossiblePairs.size() == minPairCount)
                    {
                        for (SvLinkedPair otherPair : otherPossiblePairs)
                        {
                            if (!minLinkPairs.contains(otherPair))
                                minLinkPairs.add(otherPair);
                        }
                    }
                }
            }
        }

        return minLinkPairs;
    }

    private void addPairToChain(final SvLinkedPair pair, boolean isExact)
    {
        // attempt to add to existing chain
        boolean addedToChain = false;
        boolean[] pairToChain = new boolean[2];

        SvBreakend unlinkedBeFirst = null;
        SvBreakend unlinkedBeSecond = null;
        final SvLinkedPair newPair;

        if(isExact)
        {
            newPair = pair;
        }
        else
        {
            // this pair was created from the set of possibles, but needs to make use of unlinked breakends
            unlinkedBeFirst = findUnlinkingMatchingBreakend(pair.getBreakend(true));
            unlinkedBeSecond = findUnlinkingMatchingBreakend(pair.getBreakend(false));

            if(unlinkedBeFirst == null || unlinkedBeSecond == null)
            {
                LOGGER.error("new pair breakendStart({} valid={})  and breakendEnd({} valid={}) no unlinked match found",
                        pair.getBreakend(true).toString(), unlinkedBeFirst != null,
                        pair.getBreakend(true).toString(), unlinkedBeSecond != null);
                return;
            }

            newPair = new SvLinkedPair(unlinkedBeFirst.getSV(), unlinkedBeSecond.getSV(), LINK_TYPE_TI,
                    unlinkedBeFirst.usesStart(), unlinkedBeSecond.usesStart());
        }

        for(SvChain chain : mPartialChains)
        {
            for(int be = SVI_START; be <= SVI_END; ++be)
            {
                boolean isStart = isStart(be);
                final SvVarData chainSV = chain.getChainEndSV(isStart);

                if (chain.canAddLinkedPair(newPair, isStart, true))
                {
                    if (chainSV.equals(newPair.first(), true))
                    {
                        pairToChain[SVI_START] = true;
                        pairToChain[SVI_END] = false;

                        if (chainSV != newPair.first())
                        {
                            // if the chain has a specific SV, ensure it is matched by this new pair
                            newPair.replaceFirst(chainSV);
                        }
                    }
                    else
                    {
                        pairToChain[SVI_START] = false;
                        pairToChain[SVI_END] = true;

                        if (chainSV != newPair.second())
                        {
                            newPair.replaceSecond(chainSV);
                        }
                    }

                    chain.addLink(newPair, isStart);
                    addedToChain = true;

                    LOGGER.debug("adding linked pair({} {}) to existing chain({}) {}",
                            newPair.toString(), newPair.assemblyInferredStr(), chain.id(), isStart ? "start" : "end");
                    break;
                }
            }

            if(addedToChain)
                break;
        }

        if(!addedToChain)
        {
            SvChain chain = new SvChain(mPartialChains.size());
            mPartialChains.add(chain);
            chain.addLink(newPair, true);
            pairToChain[SVI_START] = true;
            pairToChain[SVI_END] = true;

            LOGGER.debug("adding linked pair({} {}) to new chain({})",
                    newPair.toString(), newPair.assemblyInferredStr(), chain.id());
        }

        registerNewLink(newPair, pairToChain);

        if(addedToChain)
        {
            // now see if any partial chains can be linked
            reconcileChains();
        }
    }

    private SvBreakend findUnlinkingMatchingBreakend(final SvBreakend breakend)
    {
        for(SvBreakend other : mUnlinkedBreakends)
        {
            if(breakend.getOrigSV() == other.getOrigSV() && breakend.usesStart() == other.usesStart())
                return other;
        }

        return null;
    }

    private void addAssemblyLinksToChains()
    {
        if(mAssemblyLinkedPairs.isEmpty())
            return;

        for(SvLinkedPair pair : mAssemblyLinkedPairs)
        {
            addPairToChain(pair, true);
        }

        if(!mPartialChains.isEmpty())
        {
            LOGGER.debug("created {} partial chains from {} assembly links", mPartialChains.size(), mAssemblyLinkedPairs.size());
        }
    }

    private void registerNewLink(final SvLinkedPair newPair, boolean[] pairToChain)
    {
        mChainedPairs.add(newPair);

        for (int be = SVI_START; be <= SVI_END; ++be)
        {
            boolean isStart = isStart(be);

            final SvBreakend breakend = newPair.getBreakend(isStart);

            mUnlinkedSVs.remove(breakend.getSV());

            mUnlinkedBreakends.remove(breakend);

            // inefficient but need to then update the possible links by checking if any breakends for this SV
            // are still unlinked, allowing a possible link to be used
            boolean hasUnlinkedBreakend = false;

            if(mCluster.hasReplicatedSVs())
            {
                for (final SvBreakend otherBe : mUnlinkedBreakends)
                {
                    if (otherBe.getOrigSV() == breakend.getOrigSV() && otherBe.usesStart() == breakend.usesStart())
                    {
                        hasUnlinkedBreakend = true;
                        break;
                    }
                }
            }

            final SvBreakend origBreakend = breakend.getOrigSV().getBreakend(breakend.usesStart());
            List<SvLinkedPair> possibleLinks = !mSvBreakendPossibleLinks.isEmpty() ? mSvBreakendPossibleLinks.get(origBreakend) : null;

            if(possibleLinks != null)
            {
                if (!hasUnlinkedBreakend)
                {
                    removePossibleLinks(possibleLinks, breakend);
                }
            }

            // check for an opposite pairing between these 2 SVs - need to look into other breakends' lists
            final SvBreakend otherOrigBreakend = breakend.getOrigSV().getBreakend(!breakend.usesStart());

            possibleLinks = !mSvBreakendPossibleLinks.isEmpty() ? mSvBreakendPossibleLinks.get(otherOrigBreakend) : null;

            if(possibleLinks != null)
            {
                final SvBreakend otherOrigBreakendAlt = newPair.first() ==  breakend.getSV() ?
                        newPair.second().getOrigSV().getBreakend(!newPair.secondLinkOnStart()) :
                        newPair.first().getOrigSV().getBreakend(!newPair.firstLinkOnStart());

                for(SvLinkedPair pair : possibleLinks)
                {
                    if(pair.hasBreakend(otherOrigBreakend) && pair.hasBreakend(otherOrigBreakendAlt))
                    {
                        possibleLinks.remove(pair);
                        break;
                    }
                }
            }

            if(mCluster.hasReplicatedSVs())
            {
                // reduce replication counts for breakends which are added to a chain
                if (pairToChain[be])
                {
                    Integer replicationCount = mSvReplicationMap.get(breakend.getOrigSV());

                    if (replicationCount != null)
                    {
                        if (replicationCount <= 1)
                            mSvReplicationMap.remove(breakend.getOrigSV());
                        else
                            mSvReplicationMap.put(breakend.getOrigSV(), replicationCount - 1);
                    }
                }
            }
        }
    }

    private void removePossibleLinks(List<SvLinkedPair> possibleLinks, SvBreakend fullyLinkedBreakend)
    {
        if(possibleLinks == null || possibleLinks.isEmpty())
            return;

        final SvVarData linkedSV = fullyLinkedBreakend.getOrigSV();
        final SvBreakend origBreakend = fullyLinkedBreakend.getOrigSV().getBreakend(fullyLinkedBreakend.usesStart());

        int index = 0;
        while (index < possibleLinks.size())
        {
            SvLinkedPair possibleLink = possibleLinks.get(index);
            if (possibleLink.hasBreakend(linkedSV, fullyLinkedBreakend.usesStart()))
            {
                // remove this from consideration
                possibleLinks.remove(index);

                SvBreakend otherBreakend = possibleLink.getBreakend(true) == origBreakend ?
                        possibleLink.getBreakend(false) : possibleLink.getBreakend(true);

                // and remove the pair which was cached in the other breakend's possibles list
                List<SvLinkedPair> otherPossibles = mSvBreakendPossibleLinks.get(otherBreakend);

                if(otherPossibles != null)
                {
                    for (SvLinkedPair otherPair : otherPossibles)
                    {
                        if (otherPair == possibleLink)
                        {
                            otherPossibles.remove(otherPair);

                            if (otherPossibles.isEmpty())
                            {
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

        if(possibleLinks.isEmpty())
        {
            mSvBreakendPossibleLinks.remove(origBreakend);
        }
    }

    private void determinePossibleLinks()
    {
        mSvBreakendPossibleLinks.clear();

        final Map<String,List<SvBreakend>> chrBreakendMap = mCluster.getChrBreakendMap();

        for (final Map.Entry<String, List<SvBreakend>> entry : chrBreakendMap.entrySet())
        {
            final List<SvBreakend> breakendList = entry.getValue();

            for (int i = 0; i < breakendList.size() -1; ++i)
            {
                final SvBreakend lowerBreakend = breakendList.get(i);

                if(lowerBreakend.orientation() != -1)
                    continue;

                if(!isUnlinkedBreakend(lowerBreakend))
                    continue;

                for (int j = i+1; j < breakendList.size(); ++j)
                {
                    final SvBreakend upperBreakend = breakendList.get(j);

                    if(upperBreakend.orientation() != 1)
                        continue;

                    if(!isUnlinkedBreakend(upperBreakend))
                        continue;

                    if(abs(upperBreakend.position() - lowerBreakend.position()) < MIN_TEMPLATED_INSERTION_LENGTH)
                        continue;

                    // record the possible link
                    final SvVarData lowerSV = lowerBreakend.getOrigSV();
                    final SvVarData upperSV = upperBreakend.getOrigSV();

                    SvLinkedPair newPair = new SvLinkedPair(lowerSV, upperSV, LINK_TYPE_TI,
                            lowerBreakend.usesStart(), upperBreakend.usesStart());

                    List<SvLinkedPair> lowerPairs = mSvBreakendPossibleLinks.get(lowerBreakend);

                    if(lowerPairs == null)
                    {
                        lowerPairs = Lists.newArrayList();
                        mSvBreakendPossibleLinks.put(lowerBreakend, lowerPairs);
                    }

                    lowerPairs.add(newPair);

                    List<SvLinkedPair> upperPairs = mSvBreakendPossibleLinks.get(upperBreakend);

                    if(upperPairs == null)
                    {
                        upperPairs = Lists.newArrayList();
                        mSvBreakendPossibleLinks.put(upperBreakend, upperPairs);
                    }

                    upperPairs.add(newPair);
                }
            }
        }
    }

    private boolean isUnlinkedBreakend(final SvBreakend breakend)
    {
        return mUnlinkedBreakends.contains(breakend);
    }

    private void setSvReplicationCounts()
    {
        mSvReplicationMap.clear();

        if(!mCluster.hasReplicatedSVs())
            return;

        for(final SvVarData var : mCluster.getSVs())
        {
            if(var.isReplicatedSv())
                continue;

            if(var.getReplicatedCount() > 0)
            {
                mSvReplicationMap.put(var, var.getReplicatedCount());
            }
        }
    }

    private void setUnlinkedBreakends()
    {
        // make a cache of all unchained breakends in those of replicated SVs
        for(SvVarData var : mCluster.getSVs())
        {
            for (int be = SVI_START; be <= SVI_END; ++be)
            {
                boolean isStart = isStart(be);

                if (var.isNullBreakend() && !isStart)
                    continue;

                mUnlinkedBreakends.add(var.getBreakend(isStart));
            }
        }

        mUnlinkedSVs.addAll(mCluster.getSVs());
    }

    private List<SvVarData> getMaxReplicationSvIds()
    {
        List<SvVarData> maxRepIds = Lists.newArrayList();
        int maxRepCount = 1;

        for(Map.Entry<SvVarData,Integer> entry : mSvReplicationMap.entrySet())
        {
            int repCount = entry.getValue();

            if(repCount <= 1)
                continue;

            if(repCount > maxRepCount)
            {
                maxRepCount = repCount;
                maxRepIds.clear();
                maxRepIds.add(entry.getKey());
            }
            else if(repCount == maxRepCount)
            {
                if(!maxRepIds.contains(entry.getKey()))
                    maxRepIds.add(entry.getKey());
            }
        }

        return maxRepIds;
    }

    private void reconcileChains()
    {
        int index1 = 0;
        while(index1 < mPartialChains.size())
        {
            SvChain chain1 = mPartialChains.get(index1);

            boolean chainsMerged = false;

            for (int index2 = index1 + 1; index2 < mPartialChains.size(); ++index2)
            {
                SvChain chain2 = mPartialChains.get(index2);

                for (int be1 = SVI_START; be1 <= SVI_END; ++be1)
                {
                    boolean c1Start = isStart(be1);

                    for (int be2 = SVI_START; be2 <= SVI_END; ++be2)
                    {
                        boolean c2Start = isStart(be2);

                        if (chain1.canAddLinkedPair(chain2.getLinkedPair(c2Start), c1Start, false))
                        {
                            LOGGER.debug("merging chain({} links={}) with chain({} links={})",
                                    chain1.id(), chain1.getLinkCount(), chain2.id(), chain2.getLinkCount());

                            // merge chains and remove the latter
                            for (SvLinkedPair linkedPair : chain2.getLinkedPairs())
                            {
                                chain1.addLink(linkedPair, c1Start);
                            }

                            mPartialChains.remove(index2);

                            chainsMerged = true;
                            break;
                        }

                    }

                    if (chainsMerged)
                        break;
                }

                if (chainsMerged)
                    break;
            }

            if (!chainsMerged)
            {
                ++index1;
            }
        }
    }

    // OLD CHAINING METHODS
    private void findSvChainsIncrementally(final List<SvVarData> svList, List<SvChain> chainsList)
    {
        // routine flow:
        // start with any existing assembly links - these are always given priority when starting or extending chains
        // then take the open ends of the current chain and find the shortest TI from amongst the remaining unconnected SVs

        // restrictions:
        // - SVs cannot be in more than 1 chain
        // - for replicated SVs, there cannot be conflicting sets of breakend pairs (eg A-B and A-C)

        // isSpecificCluster(mCluster);
        List<SvLinkedPair> chainedPairs = Lists.newArrayList();
        List<SvLinkedPair> remainingStartLinks = Lists.newArrayList();
        remainingStartLinks.addAll(mAssemblyLinkedPairs);
        remainingStartLinks.addAll(mCluster.getInferredLinkedPairs());

        final List<SvVarData> unlinkedSvList = Lists.newArrayList();

        if(mCluster.hasSubClusters())
            unlinkedSvList.addAll(mCluster.getUnlinkedSVs());
        else
            unlinkedSvList.addAll(svList);

        List<SvLinkedPair> requiredLinks = Lists.newArrayList();

        boolean chainComplete = false;
        boolean chainLinkAdded = false;
        SvChain currentChain = null;

        while(!chainComplete || chainLinkAdded)
        {
            // start with a single linked pair
            // for each of its ends (where the first BE is labelled 'first', and the second labelled 'last'),
            // search for the closest possible linking BE from another linked pair
            // for BEs to link they must be facing (like a TI)
            if(!chainLinkAdded || currentChain == null)
            {
                // take either the next known pair or look for a new candidate
                SvLinkedPair linkedPair = null;

                if(!remainingStartLinks.isEmpty())
                {
                    linkedPair = remainingStartLinks.get(0);
                    remainingStartLinks.remove(0);
                }
                else
                {
                    linkedPair = findNewLinkedPair(chainedPairs, unlinkedSvList);

                    // nothing left to start a chain with
                    if(linkedPair == null)
                        break;
                }

                currentChain = new SvChain(chainsList.size());

                if (mLogVerbose)
                {
                    LOGGER.debug("cluster({}) starting chain({}) with linked pair({})",
                            mCluster.id(), currentChain.id(), linkedPair.toString());
                }

                currentChain.addLink(linkedPair, true);
                chainedPairs.add(linkedPair);
                reduceRemainingLists(linkedPair, unlinkedSvList, remainingStartLinks);
            }

            chainLinkAdded = false;

            // look amongst the remaining assembly links for a link to the current chain
            int index = 0;
            while(index < remainingStartLinks.size())
            {
                SvLinkedPair assemblyLink = remainingStartLinks.get(index);

                if(assemblyLink.isInferred())
                    break;

                boolean addToStart = false;
                if (currentChain.canAddLinkedPairToStart(assemblyLink))
                {
                    addToStart = true;
                    currentChain.addLink(assemblyLink, true);
                }
                else if (currentChain.canAddLinkedPairToEnd(assemblyLink))
                {
                    addToStart = false;
                    currentChain.addLink(assemblyLink, false);
                }
                else
                {
                    ++index;
                    continue;
                }

                LOGGER.debug("adding assembly linked pair({}) to chain({}) {}",
                        assemblyLink.toString(), currentChain.id(), addToStart ? "start" : "end");

                remainingStartLinks.remove(index);
                chainedPairs.add(assemblyLink);
                reduceRemainingLists(assemblyLink, unlinkedSvList, remainingStartLinks);
                chainLinkAdded = true;
            }

            // now find the closest TI from amongst the remaining unlinked variant breakends
            SvVarData chainFirstSV = currentChain.getFirstSV();
            boolean chainFirstUnlinkedOnStart = currentChain.firstLinkOpenOnStart();
            SvVarData chainLastSV = currentChain.getLastSV();
            boolean chainLastUnlinkedOnStart = currentChain.lastLinkOpenOnStart();

            requiredLinks.clear();
            requiredLinks.addAll(chainedPairs);
            requiredLinks.addAll(mAssemblyLinkedPairs); // these will restrict potential inferred TIs even if not part of a chain yet

            SvLinkedPair closestStartPair = findNextLinkedPair(requiredLinks, unlinkedSvList, chainFirstSV, chainFirstUnlinkedOnStart, 0);
            SvLinkedPair closestLastPair = findNextLinkedPair(requiredLinks, unlinkedSvList, chainLastSV, chainLastUnlinkedOnStart, 0);

            if(closestStartPair != null || closestLastPair != null)
            {
                chainLinkAdded = true;

                // don't add the same variant to both start and end - instead choose the shortest of the 2
                if (closestStartPair != null && closestLastPair != null && closestStartPair.hasAnySameVariant(closestLastPair))
                {
                    if (closestStartPair.length() < closestLastPair.length())
                        closestLastPair = null;
                    else
                        closestStartPair = null;
                }

                if (closestStartPair != null)
                {
                    if (mLogVerbose)
                    {
                        LOGGER.debug("adding linked pair({}) to chain({}) start",
                                closestStartPair.toString(), currentChain.id());
                    }

                    // add this to the chain at the start
                    currentChain.addLink(closestStartPair, true);
                    chainedPairs.add(closestStartPair);
                    reduceRemainingLists(closestStartPair, unlinkedSvList, remainingStartLinks);
                }

                if (closestLastPair != null & closestLastPair != closestStartPair)
                {
                    if (mLogVerbose)
                    {
                        LOGGER.debug("adding linked pair({}) to chain({}) end",
                                closestLastPair.toString(), currentChain.id());
                    }

                    // add this to the chain at the start
                    currentChain.addLink(closestLastPair, false);
                    chainedPairs.add(closestLastPair);
                    reduceRemainingLists(closestLastPair, unlinkedSvList, remainingStartLinks);
                }
            }

            chainComplete = (currentChain.getSvCount() == svList.size());

            if((!chainLinkAdded && closestStartPair == null && closestLastPair == null)
            || chainComplete || (remainingStartLinks.isEmpty() && unlinkedSvList.isEmpty()))
            {
                // previously: take any chain with 2 or more links or a chain from assembly links
                // if(currentChain.getLinkCount() > 1 || !currentChain.getLinkedPairs().get(0).isInferred())
                chainsList.add(currentChain);

                if(mLogVerbose)
                {
                    LOGGER.debug("cluster({}) adding {} chain({}) with {} linked pairs:",
                            mCluster.id(), chainComplete ? "complete" : "partial", currentChain.id(), currentChain.getLinkCount());
                    currentChain.logLinks();
                }

                // remove SVs (and their replicated instances) used in this chain from the unlinked list to prevent them being used in any other
                for(final SvVarData var : currentChain.getSvList())
                {
                    int i = 0;
                    while(i < unlinkedSvList.size())
                    {
                        SvVarData unlinkedSV = unlinkedSvList.get(i);
                        if(unlinkedSV.equals(var, true))
                            unlinkedSvList.remove(i);
                        else
                            ++i;
                    }
                }

                for(SvLinkedPair pair : currentChain.getLinkedPairs())
                {
                    int i = 0;
                    while(i < remainingStartLinks.size())
                    {
                        final SvLinkedPair other = remainingStartLinks.get(i);
                        if(pair.matches(other))
                            remainingStartLinks.remove(i);
                        else
                            ++i;
                    }
                }

                currentChain = null;

                if(chainComplete)
                    break;
            }
        }

        if(!chainComplete && !chainsList.isEmpty() && unlinkedSvList.size() == 1 && unlinkedSvList.get(0).type() == SGL)
        {
            tryUnlinkedSingle(chainsList, unlinkedSvList.get(0));
        }

        if (mLogVerbose)
        {
            LOGGER.debug("cluster({}) chaining finished: chains({}) unlinkedSVs({})", mCluster.id(), chainsList.size(), unlinkedSvList.size());
        }

    }

    // private static int MIN_INFERRED_CHAIN_TI_LENGTH = 100;
    private static int MIN_INFERRED_CHAIN_TI_LENGTH = MIN_TEMPLATED_INSERTION_LENGTH;

    private SvLinkedPair findNewLinkedPair(final List<SvLinkedPair> chainedPairs, final List<SvVarData> unlinkedSVs)
    {
        // try to find an inferred linked pair from which to begin a new chain
        SvLinkedPair shortestPair = null;

        for(int i = 0; i < unlinkedSVs.size(); ++i)
        {
            final SvVarData var = unlinkedSVs.get(i);

            for (int be = SVI_START; be <= SVI_END; ++be)
            {
                boolean useStart = isStart(be);

                SvLinkedPair pair = findNextLinkedPair(chainedPairs, unlinkedSVs, var, useStart, i+1);

                if(pair == null)
                    continue;

                if(shortestPair == null || pair.length() < shortestPair.length())
                    shortestPair = pair;
            }
        }

        return shortestPair;
    }

    private SvLinkedPair findNextLinkedPair(final List<SvLinkedPair> chainedPairs, final List<SvVarData> svList,
            final SvVarData var, boolean useStart, int startSvIndex)
    {
        // find the shortest templated insertion which links to this variant
        if(var.isAssemblyMatched(useStart)) // this breakend will be added when the assembly link is added
            return null;

        SvLinkedPair newPair = null;

        for(int i = startSvIndex; i < svList.size(); ++i)
        {
            final SvVarData otherVar = svList.get(i);

            if(var.equals(otherVar, true))
                continue;

            for(int be = SVI_START; be <= SVI_END; ++be)
            {
                boolean otherVarStart = isStart(be);

                if(otherVar.isNullBreakend() && !otherVarStart)
                    continue;

                if (otherVar.isAssemblyMatched(otherVarStart))
                    continue;

                if (!areLinkedSection(var, otherVar, useStart, otherVarStart, !mCluster.hasReplicatedSVs()))
                    continue;

                int tiLength = getProximity(var, otherVar, useStart, otherVarStart);

                if (tiLength <= MIN_INFERRED_CHAIN_TI_LENGTH || (newPair != null && tiLength >= newPair.length()))
                    continue;

                // check if already in a linked pair
                boolean isValid = true;
                for(final SvLinkedPair pair : chainedPairs)
                {
                    if(pair.hasBreakend(otherVar, otherVarStart))
                    {
                        isValid = false;
                        break;
                    }
                }

                if(!isValid)
                    continue;

                // form a new TI from these 2 BEs
                SvLinkedPair testPair = new SvLinkedPair(var, otherVar, SvLinkedPair.LINK_TYPE_TI, useStart, otherVarStart);

                if(mCluster.hasReplicatedSVs())
                {
                    // check that if this link has the same SVs as another link due to replicated SVs, that
                    // the breakends used in the link don't match on opposite breakends eg start-end & end-start are not allowed
                    for(final SvLinkedPair pair : chainedPairs)
                    {
                        if(!testPair.sameVariants(pair))
                            continue;

                        // if(!testPair.matches(pair) && !testPair.oppositeMatch(pair)) // previous logic was incorrect
                        if(testPair.oppositeMatch(pair))
                        {
                            isValid = false;
                            break;
                        }
                    }
                }

                if(!isValid)
                    continue;

                newPair = testPair;
            }
        }

        return newPair;
    }

    private void tryUnlinkedSingle(List<SvChain> chains, SvVarData singleSV)
    {
        // test the spare SGL against any inconsistent SVs with matching copy number change or another open SGL on end of a chain
        SvChain amendedChain = null;
        boolean addToStart = true;
        SvLinkedPair newPair = null;

        for(SvChain chain : chains)
        {
            newPair = findInconsistentSvAndSingleLink(chain, singleSV, true);

            if(newPair != null)
            {
                amendedChain = chain;
                addToStart = true;
                break;
            }

            newPair = findInconsistentSvAndSingleLink(chain, singleSV, false);

            if(newPair != null)
            {
                amendedChain = chain;
                addToStart = false;
                break;
            }

            if (chain.getFirstSV().type() == SGL && !chain.getFirstSV().equals(singleSV, true))
            {
                newPair = new SvLinkedPair(chain.getFirstSV(), singleSV, LINK_TYPE_TI, false, false);
                LOGGER.debug("adding linked pair({}) with SGL to chain({}) start", newPair.toString(), chain.id());
                amendedChain = chain;
                addToStart = true;
                break;
            }

            if (chain.getLastSV().type() == SGL && !chain.getLastSV().equals(singleSV, true))
            {
                newPair = new SvLinkedPair(chain.getLastSV(), singleSV, LINK_TYPE_TI, false, false);
                LOGGER.debug("adding linked pair({}) with SGL to chain({}) end", newPair.toString(), chain.id());
                amendedChain = chain;
                addToStart = false;
                break;
            }
        }

        if(amendedChain != null && newPair != null)
        {
            amendedChain.addLink(newPair, addToStart);

            // check if any chains can now be joined
            for(int index = 0; index < chains.size(); ++index)
            {
                final SvChain chain = chains.get(index);

                if(chain == amendedChain)
                    continue;

                if (reconcileChains(amendedChain, chain))
                {
                    LOGGER.debug("adding existing chain({}) to current chain({})", chain.id(), amendedChain.id());
                    chains.remove(index);
                    break;
                }
            }
        }
    }

    private final SvLinkedPair findInconsistentSvAndSingleLink(final SvChain chain, final SvVarData soloVar, boolean checkChainStart)
    {
        // look for a single variant with inconsistent CN change in this cluster
        // for potential pairing with a solo single
        SvVarData inconsistentVar = null;
        boolean inconsistentVarOpenOnStart = false;

        if(checkChainStart && chain.getFirstSV().hasInconsistentCopyNumberChange(chain.firstLinkOpenOnStart()))
        {
            inconsistentVar = chain.getFirstSV();
            inconsistentVarOpenOnStart = chain.firstLinkOpenOnStart();
        }
        else if(!checkChainStart && chain.getLastSV().hasInconsistentCopyNumberChange(chain.lastLinkOpenOnStart()))
        {
            inconsistentVar = chain.getLastSV();
            inconsistentVarOpenOnStart = chain.lastLinkOpenOnStart();
        }

        if(inconsistentVar == null)
            return null;

        double cnInconsistency = inconsistentVar.getSvData().ploidy() - inconsistentVar.copyNumberChange(inconsistentVarOpenOnStart);

        if(round(cnInconsistency) != round(soloVar.copyNumberChange(true)))
            return null;

        LOGGER.debug(String.format("cluster(%s) chain(%d) inconsistentSv(%s cnChg=%.2f ploidy=%.2f) joined to solo-single(%s chChg%.2f)",
                mCluster.id(), chain.id(),
                inconsistentVar.id(), inconsistentVar.getSvData().ploidy(), inconsistentVar.copyNumberChange(inconsistentVarOpenOnStart),
                soloVar.id(), soloVar.copyNumberChange(true)));

        return new SvLinkedPair(inconsistentVar, soloVar, LINK_TYPE_TI, inconsistentVarOpenOnStart, false);
    }

    private static void reduceRemainingLists(final SvLinkedPair newPair, List<SvVarData> unlinkedSvList, List<SvLinkedPair> remainingPairs)
    {
        unlinkedSvList.remove(newPair.first());
        unlinkedSvList.remove(newPair.second());

        // also remove any inferred TI which has one of these SVs
        int index = 0;
        while(index < remainingPairs.size())
        {
            final SvLinkedPair pair = remainingPairs.get(index);

            if(pair.isInferred() && pair.hasAnySameVariant(newPair))
            {
                remainingPairs.remove(index);
            }
            else
            {
                ++index;
            }
        }
    }

    private static boolean reconcileChains(SvChain firstChain, SvChain nextChain)
    {
        boolean canAddToStart = firstChain.canAddLinkedPairToStart(nextChain.getFirstLinkedPair());
        boolean canAddToEnd = firstChain.canAddLinkedPairToEnd(nextChain.getFirstLinkedPair());

        if(canAddToStart || canAddToEnd)
        {
            for(SvLinkedPair linkedPair : nextChain.getLinkedPairs())
            {
                firstChain.addLink(linkedPair, canAddToStart);
            }

            return true;
        }

        canAddToStart = firstChain.canAddLinkedPairToStart(nextChain.getLastLinkedPair());
        canAddToEnd = firstChain.canAddLinkedPairToEnd(nextChain.getLastLinkedPair());

        if(canAddToStart || canAddToEnd)
        {
            // add in reverse
            for(int index = nextChain.getLinkedPairs().size() - 1; index >= 0; --index)
            {
                SvLinkedPair linkedPair = nextChain.getLinkedPairs().get(index);
                firstChain.addLink(linkedPair, canAddToStart);
            }

            return true;
        }

        return false;
    }

}
