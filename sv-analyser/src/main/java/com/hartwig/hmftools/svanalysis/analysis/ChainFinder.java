package com.hartwig.hmftools.svanalysis.analysis;

import static java.lang.Math.abs;
import static java.lang.Math.round;

import static com.hartwig.hmftools.common.variant.structural.StructuralVariantType.SGL;
import static com.hartwig.hmftools.svanalysis.analysis.LinkFinder.MIN_TEMPLATED_INSERTION_LENGTH;
import static com.hartwig.hmftools.svanalysis.analysis.LinkFinder.areLinkedSection;
import static com.hartwig.hmftools.svanalysis.analysis.SvUtilities.getProximity;
import static com.hartwig.hmftools.svanalysis.types.SvChain.checkChainReplication;
import static com.hartwig.hmftools.svanalysis.types.SvCluster.isSpecificCluster;
import static com.hartwig.hmftools.svanalysis.types.SvVarData.SVI_END;
import static com.hartwig.hmftools.svanalysis.types.SvVarData.SVI_START;
import static com.hartwig.hmftools.svanalysis.types.SvVarData.isStart;
import static com.hartwig.hmftools.svanalysis.types.SvLinkedPair.LINK_TYPE_TI;

import java.util.HashMap;
import java.util.List;
import java.util.Map;

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

    private List<SvChain> mPartialChains;
    private int mNextChainId;
    private List<SvVarData> mFullyLinkedSVs;
    private List<SvVarData> mUnlinkedSVs;
    private Map<SvBreakend,List<SvBreakend>> mUnlinkedBreakendMap;
    private Map<SvVarData,Integer> mSvReplicationMap;

    private Map<SvBreakend,List<SvLinkedPair>> mSvBreakendPossibleLinks;

    private boolean mIsValid;
    private boolean mLogVerbose;
    private boolean mLogWorking;

    public ChainFinder()
    {
        mNextChainId = 0;
        mPartialChains = Lists.newArrayList();
        mUnlinkedSVs = Lists.newArrayList();
        mUnlinkedBreakendMap = new HashMap();
        mFullyLinkedSVs = Lists.newArrayList();
        mSvReplicationMap = new HashMap();
        mSvBreakendPossibleLinks = new HashMap();

        mLogVerbose = false;
        mLogWorking = false;
        mIsValid = true;
    }

    public void initialise(SvCluster cluster)
    {
        mNextChainId = 0;
        mCluster = cluster;
        mIsValid = true;

        mAssemblyLinkedPairs = Lists.newArrayList(cluster.getAssemblyLinkedPairs());

        mPartialChains.clear();
        mFullyLinkedSVs.clear();
        mUnlinkedSVs.clear();
        mUnlinkedBreakendMap.clear();
    }

    public void setLogVerbose(boolean toggle) { mLogVerbose = toggle; }
    public void setLogWorking(boolean toggle) { mLogWorking = toggle; }

    public void formClusterChains(boolean assembledLinksOnly)
    {
        List<SvVarData> svList = Lists.newArrayList(mCluster.getSVs());

        if(svList.size() < 2)
            return;

        if (mCluster.getCount() >= 4)
        {
            LOGGER.debug("cluster({}) assemblyLinks({}) svCount({} rep={})",
                    mCluster.id(), mAssemblyLinkedPairs.size(), mCluster.getUniqueSvCount(), mCluster.getCount());
        }

        isSpecificCluster(mCluster);
        // mLogWorking = isSpecificCluster(mCluster);

        buildSvChains(assembledLinksOnly);

        if(!mIsValid)
        {
            LOGGER.warn("cluster({}) chain finding failed", mCluster.id());
            return;
        }

        mCluster.getChains().clear();

        // add these chains to the cluster, but skip any which are identical to existing ones,
        // which can happen for clusters with replicated SVs
        mPartialChains.stream().forEach(chain -> checkAddNewChain(chain));

        for(int i = 0; i < mCluster.getChains().size(); ++i)
        {
            final SvChain chain = mCluster.getChains().get(i);
            LOGGER.debug("cluster({}) added chain({}) with {} linked pairs:",
                    mCluster.id(), chain.id(), chain.getLinkCount());
            chain.logLinks();
            chain.setId(i); // set after logging so can compare with logging during building
        }
    }

    private void checkAddNewChain(final SvChain newChain)
    {
        if(!mCluster.hasReplicatedSVs() || mCluster.getChains().isEmpty())
        {
            mCluster.addChain(newChain, false);
            return;
        }

        // any identical chains will have their replicated SVs entirely removed
        for(final SvChain chain : mCluster.getChains())
        {
            if(chain.identicalChain(newChain))
            {
                boolean allReplicatedSVs = newChain.getUniqueSvCount() == 0;

                LOGGER.debug("cluster({}) skipping duplicate chain({}) vs origChain({}) all replicated({})",
                        mCluster.id(), newChain.id(), chain.id(), allReplicatedSVs);

                // remove these replicated SVs as well as the replicated chain
                if(allReplicatedSVs)
                {
                    for (final SvVarData var : newChain.getSvList())
                    {
                        mCluster.removeReplicatedSv(var);
                    }
                }

                return;
            }
        }

        mCluster.addChain(newChain, false);
    }

    private void buildSvChains(boolean assembledLinksOnly)
    {
        setUnlinkedBreakends();

        // first make chains out of any assembly links
        addAssemblyLinksToChains();

        if(!assembledLinksOnly)
        {
            setSvReplicationCounts();

            determinePossibleLinks();

            // now find the best next candidate link giving priority to replication count, following by min link resolution
            // and lastly shortest distance

            while (true)
            {
                // first check if there are SVs with a higher replication count, and if so favour these first
                List<SvVarData> maxRepSVs = !mSvReplicationMap.isEmpty() ? getMaxReplicationSvIds() : null;

                List<SvBreakend> breakendList = Lists.newArrayList();

                for (final SvBreakend breakend : mUnlinkedBreakendMap.keySet())
                {
                    if (maxRepSVs != null && !maxRepSVs.contains(breakend.getSV()))
                        continue;

                    breakendList.add(breakend);
                }

                boolean isMaxReplicated = maxRepSVs != null && !maxRepSVs.isEmpty();

                if (mLogWorking && isMaxReplicated)
                {
                    for (SvVarData var : maxRepSVs)
                    {
                        LOGGER.debug("restricted to rep SV: {} repCount({})", var.id(), mSvReplicationMap.get(var));
                    }
                }

                // next take the pairings with the least alternatives
                List<SvLinkedPair> possiblePairs = findRestrictedPairs(breakendList, isMaxReplicated);

                if (possiblePairs.isEmpty())
                {
                    if (isMaxReplicated)
                    {
                        // these high-replication SVs yielded no possible links so remove them from consideration
                        for (final SvVarData var : maxRepSVs)
                        {
                            if (mLogVerbose)
                            {
                                LOGGER.debug("cluster({}) removing high-replicated SV({} {})", mCluster.id(), var.posId(), var.type());
                            }

                            mSvReplicationMap.remove(var);
                        }

                        continue;
                    }

                    break;
                }

                processPossiblePairs(possiblePairs, isMaxReplicated);
                checkProgress();
            }
        }

        if(mLogVerbose)
        {
            int breakendCount = (int)mUnlinkedBreakendMap.values().stream().count();

            List<SvVarData> uniqueUnlinkedSVs = Lists.newArrayList();

            for(final SvVarData var : mUnlinkedSVs)
            {
                if(!uniqueUnlinkedSVs.contains(var.getOrigSV()))
                    uniqueUnlinkedSVs.add(var.getOrigSV());
            }

            LOGGER.debug("cluster({}) chaining finished: chains({}) unlinked SVs({} unique={}) breakends({} reps={})",
                    mCluster.id(), mPartialChains.size(), mUnlinkedSVs.size(), uniqueUnlinkedSVs.size(),
                    mUnlinkedBreakendMap.size(), breakendCount);
        }
    }

    private void processPossiblePairs(List<SvLinkedPair> possiblePairs, boolean isMaxReplicated)
    {
        while (!possiblePairs.isEmpty())
        {
            SvLinkedPair shortestPair = null;
            for (SvLinkedPair pair : possiblePairs)
            {
                if (mLogWorking)
                {
                    LOGGER.debug("possible pair: {} length({})", pair.toString(), pair.length());
                }

                if (shortestPair == null || pair.length() < shortestPair.length())
                {
                    shortestPair = pair;
                }
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

            if (mLogWorking)
            {
                LOGGER.debug("shortest possible pair: {} length({})", shortestPair.toString(), shortestPair.length());
            }

            int pairRepeatCount = 1;

            if(isMaxReplicated)
            {
                List<SvBreakend> beListStart = mUnlinkedBreakendMap.get(shortestPair.getBreakend(true));
                List<SvBreakend> beListEnd = mUnlinkedBreakendMap.get(shortestPair.getBreakend(false));

                if(beListStart != null && beListStart.size() > 1 && beListEnd != null && beListStart.size() == beListEnd.size())
                {
                    pairRepeatCount = beListStart.size();

                    LOGGER.debug("cluster({}) repeating pair({}) {} times",
                            mCluster.id(), shortestPair.toString(), pairRepeatCount);
                }
            }

            for(int i = 0; i < pairRepeatCount; ++i)
            {
                // the actual pair being added by be drawn from unlinked breakends, not the real ones
                addPairToChain(shortestPair, false);
            }

            if(!mIsValid)
                return;
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
        boolean[] pairToChain = {false, false};

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
            unlinkedBeFirst = findUnlinkedMatchingBreakend(pair.getBreakend(true));
            unlinkedBeSecond = findUnlinkedMatchingBreakend(pair.getBreakend(false));

            if(unlinkedBeFirst == null || unlinkedBeSecond == null)
            {
                mIsValid = false;
                LOGGER.error("new pair breakendStart({} valid={}) and breakendEnd({} valid={}) no unlinked match found",
                        pair.getBreakend(true).toString(), unlinkedBeFirst != null,
                        pair.getBreakend(false).toString(), unlinkedBeSecond != null);
                return;
            }

            newPair = new SvLinkedPair(unlinkedBeFirst.getSV(), unlinkedBeSecond.getSV(), LINK_TYPE_TI,
                    unlinkedBeFirst.usesStart(), unlinkedBeSecond.usesStart());
        }

        for(SvChain chain : mPartialChains)
        {
            // test this link against both ends to a chain since closed chains are not allowed
            boolean addToStart = false;
            for(int be = SVI_START; be <= SVI_END; ++be)
            {
                boolean isStart = isStart(be);
                final SvVarData chainSV = chain.getChainEndSV(isStart);

                if (chain.canAddLinkedPair(newPair, isStart, true))
                {
                    addToStart = isStart;

                    if (chainSV.equals(newPair.first(), true))
                    {
                        pairToChain[SVI_START] = true;
                    }
                    else
                    {
                        pairToChain[SVI_END] = true;
                    }
                }
            }

            if(pairToChain[SVI_START] && pairToChain[SVI_END])
            {
                LOGGER.debug("skipping linked pair({}) would close existing chain({})", newPair.toString(), chain.id());
                addedToChain = true;
            }
            else if(!pairToChain[SVI_START] && !pairToChain[SVI_END])
            {
                continue;
            }
            else
            {
                final SvVarData chainSV = addToStart ? chain.getFirstSV() : chain.getLastSV();

                if (pairToChain[SVI_START] && chainSV != newPair.first())
                {
                    // if the chain has a specific SV, ensure it is matched by this new pair
                    newPair.replaceFirst(chainSV);
                }
                else if (pairToChain[SVI_END] && chainSV != newPair.second())
                {
                    newPair.replaceSecond(chainSV);
                }

                chain.addLink(newPair, addToStart);
                addedToChain = true;

                LOGGER.debug("adding linked pair({} {} len={}) to existing chain({}) {}",
                        newPair.toString(), newPair.assemblyInferredStr(), newPair.length(), chain.id(), addToStart ? "start" : "end");
            }

            if(addedToChain)
                break;
        }

        if(!addedToChain)
        {
            SvChain chain = new SvChain(mNextChainId++);
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

    private SvBreakend findUnlinkedMatchingBreakend(final SvBreakend breakend)
    {
        // get the next available breakend (thereby reducing the replicated instances)
        final List<SvBreakend> breakendList = mUnlinkedBreakendMap.get(breakend.getOrigBreakend());

        if(breakendList == null || breakendList.isEmpty())
            return null;

        return breakendList.get(0);
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
        for (int be = SVI_START; be <= SVI_END; ++be)
        {
            boolean isStart = isStart(be);

            final SvBreakend breakend = newPair.getBreakend(isStart);

            mUnlinkedSVs.remove(breakend.getSV());

            final SvBreakend origBreakend = breakend.getOrigBreakend();

            final List<SvBreakend> breakendList = mUnlinkedBreakendMap.get(origBreakend);

            if(breakendList == null || breakendList.isEmpty())
            {
                LOGGER.error("breakend({}) list already empty", origBreakend.toString());
                mIsValid = false;
                return;
            }

            breakendList.remove(breakend);

            boolean hasUnlinkedBreakend = true;
            if(breakendList.isEmpty())
            {
                mUnlinkedBreakendMap.remove(origBreakend);
                hasUnlinkedBreakend = false;

                if(mLogWorking)
                {
                    LOGGER.debug("breakend({}) as no more possible links", breakend);
                }
            }

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

            possibleLinks = otherOrigBreakend != null && !mSvBreakendPossibleLinks.isEmpty() ? mSvBreakendPossibleLinks.get(otherOrigBreakend) : null;

            if(possibleLinks != null)
            {
                final SvBreakend otherOrigBreakendAlt = newPair.first() ==  breakend.getSV() ?
                        newPair.second().getOrigSV().getBreakend(!newPair.secondLinkOnStart()) :
                        newPair.first().getOrigSV().getBreakend(!newPair.firstLinkOnStart());

                if(otherOrigBreakendAlt != null)
                {
                    for (SvLinkedPair pair : possibleLinks)
                    {
                        if (pair.hasBreakend(otherOrigBreakend) && pair.hasBreakend(otherOrigBreakendAlt))
                        {
                            possibleLinks.remove(pair);
                            break;
                        }
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
                        if (replicationCount <= 2)
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
        final SvBreakend origBreakend = fullyLinkedBreakend.getOrigBreakend();

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
                                if(mLogWorking)
                                {
                                    LOGGER.debug("breakend({}) as no more possible links", otherBreakend);
                                }

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
            if(mLogWorking)
            {
                LOGGER.debug("breakend({}) as no more possible links", origBreakend);
            }

            mSvBreakendPossibleLinks.remove(origBreakend);
        }
    }

    private void determinePossibleLinks()
    {
        // form a map of each breakend to its set of all other breakends which can form a valid TI
        // need to exclude breakends which are already assigned to an assembled TI
        // unless replication permits additional instances of it
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

                if(alreadyLinkedBreakend(lowerBreakend))
                    continue;

                for (int j = i+1; j < breakendList.size(); ++j)
                {
                    final SvBreakend upperBreakend = breakendList.get(j);

                    if(upperBreakend.orientation() != 1)
                        continue;

                    if(upperBreakend.getSV() == lowerBreakend.getSV())
                        continue;

                    if(alreadyLinkedBreakend(upperBreakend))
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

    private boolean alreadyLinkedBreakend(final SvBreakend breakend)
    {
        Integer beRepCount = null;
        int beRepCountRemainder = 0;

        for(SvLinkedPair pair : mAssemblyLinkedPairs)
        {
            if(pair.hasBreakend(breakend, true))
            {
                // check whether replication would still allow this breakend to be used again
                if(beRepCount == null)
                {
                    beRepCount = mSvReplicationMap.get(breakend.getOrigSV());

                    if (beRepCount == null)
                        return true;

                    beRepCountRemainder = beRepCount;
                }

                // each time this breakend is connected to another breakend via assembly,
                // it subtracts from it potential usage for inferred links
                final SvBreakend otherBreakend = pair.getOtherBreakend(breakend);
                Integer otherBeRepCount = mSvReplicationMap.get(otherBreakend.getOrigSV());

                if(otherBeRepCount == null)
                    --beRepCountRemainder;
                else
                    beRepCountRemainder -= otherBeRepCount;

                if(beRepCountRemainder <= 0)
                    return true;
            }
        }

        return false;
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
        for(final SvVarData var : mCluster.getSVs())
        {
            for (int be = SVI_START; be <= SVI_END; ++be)
            {
                boolean isStart = isStart(be);

                if (var.isNullBreakend() && !isStart)
                    continue;

                final SvBreakend breakend = var.getBreakend(isStart);
                final SvBreakend origBreakend = breakend.getOrigBreakend();

                List<SvBreakend> breakends = mUnlinkedBreakendMap.get(origBreakend);

                if(breakends == null)
                {
                    breakends = Lists.newArrayList();
                    mUnlinkedBreakendMap.put(origBreakend, breakends);
                }

                breakends.add(breakend);
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

                            if(c2Start)
                            {
                                // merge chains and remove the latter
                                for (SvLinkedPair linkedPair : chain2.getLinkedPairs())
                                {
                                    chain1.addLink(linkedPair, c1Start);
                                }
                            }
                            else
                            {
                                // add in reverse
                                for (int index = chain2.getLinkedPairs().size() - 1; index >= 0; --index)
                                {
                                    SvLinkedPair linkedPair = chain2.getLinkedPairs().get(index);
                                    chain1.addLink(linkedPair, c1Start);
                                }
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

    private void checkProgress()
    {
        if(mCluster.getCount() < 100)
            return;

        if((mPartialChains.size() % 100) == 0 || (mUnlinkedBreakendMap.size() % 100) == 0)
        {
            LOGGER.debug("cluster({}) progress: SVs({}) partialChains({}) unlinked(SVs={} breakends={}) replicatedSVs({})",
                    mCluster.id(), mCluster.getCount(), mPartialChains.size(), mUnlinkedSVs.size(),
                    mUnlinkedBreakendMap.size(), mSvReplicationMap.size());
        }
    }

}
