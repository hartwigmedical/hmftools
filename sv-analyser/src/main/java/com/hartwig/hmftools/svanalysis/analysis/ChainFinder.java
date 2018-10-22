package com.hartwig.hmftools.svanalysis.analysis;

import static java.lang.Math.abs;
import static java.lang.Math.max;
import static java.lang.Math.min;
import static java.lang.Math.round;

import static com.hartwig.hmftools.common.variant.structural.StructuralVariantType.INS;
import static com.hartwig.hmftools.common.variant.structural.StructuralVariantType.INV;
import static com.hartwig.hmftools.svanalysis.analysis.LinkFinder.areLinkedSection;
import static com.hartwig.hmftools.svanalysis.analysis.SvUtilities.getProximity;
import static com.hartwig.hmftools.svanalysis.types.SvVarData.RELATION_TYPE_NEIGHBOUR;
import static com.hartwig.hmftools.svanalysis.types.SvVarData.SVI_END;
import static com.hartwig.hmftools.svanalysis.types.SvVarData.SVI_START;
import static com.hartwig.hmftools.svanalysis.types.SvVarData.isStart;
import static com.hartwig.hmftools.svanalysis.types.SvLinkedPair.LINK_TYPE_SGL;
import static com.hartwig.hmftools.svanalysis.types.SvLinkedPair.LINK_TYPE_TI;

import java.util.List;
import java.util.Map;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.utils.PerformanceCounter;
import com.hartwig.hmftools.svanalysis.types.SvCNData;
import com.hartwig.hmftools.svanalysis.types.SvChain;
import com.hartwig.hmftools.svanalysis.types.SvCluster;
import com.hartwig.hmftools.svanalysis.types.SvVarData;
import com.hartwig.hmftools.svanalysis.types.SvLinkedPair;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;

public class ChainFinder
{
    private static double MIN_CHAIN_PERCENT = 0.25;
    private static double LOW_PLOIDY_LIMIT = 0.05;
    private static int MAX_COMPLETE_CHAINS = 20;
    private static int MAX_PATH_ITERATIONS = 250000;

    private static final Logger LOGGER = LogManager.getLogger(ChainFinder.class);

    final SvUtilities mUtils;
    private String mSampleId;
    private SvCluster mCluster;
    private List<SvLinkedPair> mAssemblyLinkedPairs;
    private List<SvLinkedPair> mInferredLinkedPairs;
    private List<SvChain> mCompleteChains;
    private List<SvChain> mIncompleteChains;
    private boolean mHasExistingChains;
    private boolean mRequireFullChains;
    private int mMinIncompleteChainCount;
    private int mReqChainCount;
    private int mPathInterations;
    private boolean mLogVerbose;
    private boolean mLogPathFinding;

    PerformanceCounter mContinuousFinderPc;
    PerformanceCounter mRecursiveFinderPc;

    public ChainFinder(final SvUtilities utils)
    {
        mUtils = utils;
        mRequireFullChains = false;
        mCompleteChains = Lists.newArrayList();
        mIncompleteChains = Lists.newArrayList();
        mLogVerbose = false;
        mPathInterations = 0;
        mReqChainCount = 0;
        mMinIncompleteChainCount = 0;
        mLogPathFinding = false;

        mContinuousFinderPc = new PerformanceCounter("Continuous");
        mRecursiveFinderPc = new PerformanceCounter("Recursive");
    }

    public final PerformanceCounter getContinuousFinderPc() { return mContinuousFinderPc; }
    public final PerformanceCounter getRecursiveFinderPc() { return mRecursiveFinderPc; }

    public void initialise(final String sampleId, SvCluster cluster)
    {
        mSampleId = sampleId;
        mCluster = cluster;
        mHasExistingChains = !mCluster.getChains().isEmpty();
        mPathInterations = 0;
        mReqChainCount = 0;
        mMinIncompleteChainCount = 0;

        mAssemblyLinkedPairs = Lists.newArrayList();

        for(SvLinkedPair pair : cluster.getAssemblyLinkedPairs())
        {
            boolean alreadyLinked = false;
            for(final SvChain chain : cluster.getChains())
            {
                if (chain.hasLinkedPair(pair))
                {
                    alreadyLinked = true;
                    break;
                }
            }

            if(!alreadyLinked)
                mAssemblyLinkedPairs.add(pair);
        }

        mInferredLinkedPairs = Lists.newArrayList();

        for(SvLinkedPair pair : cluster.getInferredLinkedPairs())
        {
            boolean alreadyLinked = false;
            for(final SvChain chain : cluster.getChains())
            {
                if (chain.hasLinkedPair(pair))
                {
                    alreadyLinked = true;
                    break;
                }
            }

            if(!alreadyLinked)
                mInferredLinkedPairs.add(pair);
        }

        mCompleteChains.clear();
        mIncompleteChains.clear();
    }

    public void setRequireFullChains(boolean toggle) { mRequireFullChains = toggle; }
    public void setLogVerbose(boolean toggle) { mLogVerbose = toggle; }

    public final List<SvChain> getCompleteChains() { return mCompleteChains; }

    public boolean formClusterChains()
    {
        // take the assembly links as a given and then try out the inferred links to see if a single chain can be formed from all the breakends

        // only factor in templated insertions to form chains, even if they also exhibit DBs
        List<SvLinkedPair> inferredTIs = collectRelevantLinkedPairs();

        List<SvVarData> svList = collectRelevantSVs();

        mReqChainCount = svList.size();

        if(!mRequireFullChains)
            mMinIncompleteChainCount = (int)round(MIN_CHAIN_PERCENT * mReqChainCount);

        // check whether there are enough potential links to form full or near-to-full chains
        if (!mCluster.hasSubClusters() && !checkLinksPotential(svList, inferredTIs))
            return false;

        if (mCluster.getCount() >= 4)
        {
            LOGGER.debug("sample({}) cluster({}) assemblyLinks({}) inferredLinks({} of total={}) svCount({} all={}) existingChains({})",
                    mSampleId, mCluster.getId(), mAssemblyLinkedPairs.size(), inferredTIs.size(), mInferredLinkedPairs.size(),
                    svList.size(), mCluster.getCount(), mCluster.getChains().size());
        }

        boolean hasChains = formChains(svList, inferredTIs);

        return hasChains;
    }

    private List<SvLinkedPair> collectRelevantLinkedPairs()
    {
        List<SvLinkedPair> inferredTIs = Lists.newArrayList();
        for (SvLinkedPair pair : mInferredLinkedPairs)
        {
            if (!(pair.linkType() == LINK_TYPE_TI || pair.linkType() == LINK_TYPE_SGL))
                continue;

            /*
            if(mRequireFullChains)
            {
                // check for similarity of copy number change
                double cnc1 = pair.first().copyNumberChange(pair.firstLinkOnStart());
                double cnc2 = pair.second().copyNumberChange(pair.secondLinkOnStart());
                double copyNumberChangeDiff = abs(cnc1 - cnc2);
                double copyNumberChangeDiffPerc = copyNumberChangeDiff / max(abs(cnc1), abs(cnc2));

                if (copyNumberChangeDiff > MAX_COPY_NUMBER_DIFF && copyNumberChangeDiffPerc > MAX_COPY_NUMBER_DIFF_PERC)
                    continue;
            }
            */

            inferredTIs.add(pair);
        }

        return inferredTIs;
    }

    private List<SvVarData> collectRelevantSVs()
    {
        List<SvVarData> svList = Lists.newArrayList();

            // remove any non-applicable or suspect SVs
        for(final SvVarData var : mCluster.getSVs())
        {
            if(var.getSvData().ploidy() < LOW_PLOIDY_LIMIT)
                continue;

            if(var.type() == INS)
                continue;

            // skip simple DELs, DUPs and INS even if in assembly links
            if(var.isSimpleType() && var.getNearestSvRelation() == RELATION_TYPE_NEIGHBOUR)
                continue;

            svList.add(var);
        }

        return svList;
    }

    private boolean checkLinksPotential(List<SvVarData> svList, List<SvLinkedPair> inferredLinkedPairs)
    {
        int minRequiredLinks = svList.size() - 1;

        if(!mRequireFullChains)
            minRequiredLinks = mMinIncompleteChainCount - 1;

        if (mAssemblyLinkedPairs.size() + inferredLinkedPairs.size() < minRequiredLinks && !mHasExistingChains)
        {
            LOGGER.debug("cluster({}) insufficient links(assembly={} inferred={}) vs svCount({}) to form chains",
                    mCluster.getId(), mAssemblyLinkedPairs.size(), inferredLinkedPairs.size(), svList.size());
            return false;
        }

        // count up unique linked-pair breakends to see whether there is the potential to form sufficiently long chains
        int singleBEMatchCount = 0;
        int bothBEMatchCount = 0;

        // factor in those SVs already chained
        for(final SvChain chain : mCluster.getChains())
        {
            bothBEMatchCount += chain.getSvCount();
        }

        List<SvLinkedPair> allLinks = Lists.newArrayList();
        allLinks.addAll(mAssemblyLinkedPairs);
        allLinks.addAll(inferredLinkedPairs);

        for (final SvVarData var : svList)
        {
            boolean startMatched = false;
            boolean endMatched = false;

            for (final SvLinkedPair linkedPair : allLinks)
            {
                if(var.isNullBreakend() && linkedPair.hasVariantBE(var, false))
                {
                    ++bothBEMatchCount;
                    continue;
                }

                if (linkedPair.hasVariantBE(var, true))
                    startMatched = true;
                else if (linkedPair.hasVariantBE(var, false))
                    endMatched = true;

                if (startMatched && endMatched)
                {
                    ++bothBEMatchCount;
                    break;
                }
            }

            if(!(startMatched && endMatched) && (startMatched || endMatched))
                ++singleBEMatchCount;
        }

        if(bothBEMatchCount >= minRequiredLinks || (singleBEMatchCount <= 2 && bothBEMatchCount >= minRequiredLinks - 1))
            return true;

        LOGGER.debug("cluster({}) insufficient matchedVars(both={} single={}) vs svCount({}) from links(assembly={} inferred={}) to form chains",
                mCluster.getId(), bothBEMatchCount, singleBEMatchCount, svList.size(), mAssemblyLinkedPairs.size(), inferredLinkedPairs.size());
        return false;
    }

    private boolean formChains(List<SvVarData> svList, List<SvLinkedPair> inferredLinkedPairs)
    {
        // find all the combinations of linked pairs where every breakend end except at most 2 are covered by a linked pair
        if(mHasExistingChains)
        {
            // partialChains.addAll(mCluster.getChains());
            // linkChains(partialChains, mReqChainCount);

            if(!mRequireFullChains)
            {
                int maxSubChainCount = 0;
                for(final SvChain chain : mCluster.getChains())
                {
                    LOGGER.debug("cluster({}) has existing chain({} svs={})", mCluster.getId(), chain.getId(), chain.getSvCount());

                    maxSubChainCount = max(maxSubChainCount, chain.getSvCount());
                }

                mMinIncompleteChainCount = max(maxSubChainCount, mMinIncompleteChainCount);
            }
        }

        mContinuousFinderPc.start();
        findContinuousChains(svList, inferredLinkedPairs);
        mContinuousFinderPc.stop();

        // first check for any complete chains, and if found, add the shortest one
        if(cacheCompleteChain())
            return true;

        // otherwise add the longest mutually exclusive chains
        cacheIncompleteChains();

        return false;
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
                continue; // only consider reciprocal chains if this short

            if (shortestFullChain == null || chain.getLength() < shortestLength)
            {
                shortestFullChain = chain;
                shortestLength = chain.getLength();
            }
        }

        if (shortestFullChain == null)
            return false;

        shortestFullChain.setId(mCluster.getChains().size());

        if(shortestFullChain.getLinkCount() >= 3)
        {
            LOGGER.info("sample({}) cluster({}) adding complete chain({}) length({}) with {} linked pairs",
                    mSampleId, mCluster.getId(), shortestFullChain.getId(), shortestFullChain.getLength(), shortestFullChain.getLinkCount());

            shortestFullChain.logLinks();
        }

        mCluster.getChains().clear();
        mCluster.addChain(shortestFullChain);
        mCluster.setIsFullyChained(true);
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

            if(maxLengthChain.getLinkCount() >= 3)
            {
                LOGGER.debug("sample({}) cluster({}) adding incomplete chain({}) length({}) with {} linked pairs",
                        mSampleId, mCluster.getId(), maxLengthChain.getId(), maxLengthChain.getLength(), maxLengthChain.getLinkCount());

                maxLengthChain.logLinks();
            }

            mCluster.addChain(maxLengthChain);
            chainedSVs.addAll(maxLengthChain.getSvList());
            mIncompleteChains.remove(maxLengthChain);
        }
    }

    public void findContinuousChains(final List<SvVarData> svList, List<SvLinkedPair> inferredLinkedPairs)
    {
        List<SvChain> chains = Lists.newArrayList();

        // re-include assembly links from any sub-clusters and all existing chains will be reformed
        if(mCluster.hasSubClusters())
        {
            for(SvCluster subCluster : mCluster.getSubClusters())
            {
                if(!subCluster.isFullyChained())
                    mAssemblyLinkedPairs.addAll(subCluster.getAssemblyLinkedPairs());

            }
        }

        findSvChainsIncrementally(svList, chains, inferredLinkedPairs);

        if(chains.size() == 1 && chains.get(0).getSvCount() == mReqChainCount)
        {
            SvChain completeChain = chains.get(0);
            mCompleteChains.add(completeChain);
            return;
        }

        for(SvChain chain : chains)
        {
            if (chain.getSvCount() < 2)
                continue;

            mIncompleteChains.add(chain);

            if(mLogVerbose)
            {
                LOGGER.debug("cluster({}) found incomplete chain({} svs={})", mCluster.getId(), chain.getId(), chain.getSvCount());
            }
        }
    }

    public void findSvChainsIncrementally(final List<SvVarData> svList, List<SvChain> chainsList, List<SvLinkedPair> inferredLinkedPairs)
    {
        // routine flow:
        // start with any existing partial chains and unconnected assembly links - these are always used
        // first to either start a new chain or link to existing chains
        // then take the open ends of the current chain and find the shortest templated insertion from amongst the remaining unconnected SVs

        List<SvChain> partialChains = Lists.newArrayList();
        partialChains.addAll(mCluster.getChains());

        if(mAssemblyLinkedPairs.isEmpty() && inferredLinkedPairs.isEmpty() && partialChains.isEmpty())
            return;

        List<SvLinkedPair> chainedPairs = Lists.newArrayList();
        List<SvLinkedPair> remainingStartLinks = Lists.newArrayList();
        remainingStartLinks.addAll(mAssemblyLinkedPairs);
        remainingStartLinks.addAll(inferredLinkedPairs);

        final List<SvVarData> unlinkedSvList = Lists.newArrayList();

        if(mCluster.hasSubClusters())
            unlinkedSvList.addAll(mCluster.getUnlinkedSVs());
        else
            unlinkedSvList.addAll(svList);

        // remove links and SVs already allocated to partial chains, but add the open ends of chains to the unlinked list
        for(final SvChain chain : partialChains)
        {
            chainedPairs.addAll(chain.getLinkedPairs());

            if(!unlinkedSvList.contains(chain.getFirstSV()))
                unlinkedSvList.add(chain.getFirstSV());

            if(!unlinkedSvList.contains(chain.getLastSV()))
                unlinkedSvList.add(chain.getLastSV());
        }

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

            if(!chainLinkAdded)
            {
                if(!partialChains.isEmpty())
                {
                    currentChain = new SvChain(partialChains.get(0));
                    partialChains.remove(0);

                    if (mLogVerbose)
                    {
                        LOGGER.debug("sample({}) cluster({}) building from existng chain({}) with {} SVs",
                                mSampleId, mCluster.getId(), currentChain.getId(), currentChain.getSvCount());
                    }
                }
                else if(!remainingStartLinks.isEmpty())
                {
                    SvLinkedPair linkedPair = remainingStartLinks.get(0);
                    remainingStartLinks.remove(0);

                    currentChain = new SvChain(chainsList.size());

                    if (mLogVerbose)
                    {
                        LOGGER.debug("sample({}) cluster({}) starting chain({}) with linked pair({})",
                                mSampleId, mCluster.getId(), currentChain.getId(), linkedPair.toString());
                    }

                    currentChain.addLink(linkedPair, true);
                    chainedPairs.add(linkedPair);
                    reduceRemainingLists(linkedPair, unlinkedSvList, remainingStartLinks);
                }
                else
                {
                    // nothing left to start a chain with
                    break;
                }
            }

            chainLinkAdded = false;

            // see if any of the remaining partial chains can now be joined
            int chainIndex = 0;
            while (chainIndex < partialChains.size())
            {
                final SvChain partialChain = partialChains.get(chainIndex);

                if (reconcileChains(currentChain, partialChain))
                {
                    LOGGER.debug("adding existing chain({}) to current chain({})", partialChain.getId(), currentChain.getId());
                    chainLinkAdded = false;
                    partialChains.remove(chainIndex);
                }
                else
                {
                    ++chainIndex;
                }
            }

            // look amongst the remaining assembly links for a link to the current chain
            int index = 0;
            while(index < remainingStartLinks.size())
            {
                SvLinkedPair assemblyLink = remainingStartLinks.get(index);

                if(assemblyLink.isInferred())
                    break;

                if (currentChain.canAddLinkedPairToStart(assemblyLink))
                {
                    currentChain.addLink(assemblyLink, true);
                }
                else if (currentChain.canAddLinkedPairToEnd(assemblyLink))
                {
                    currentChain.addLink(assemblyLink, false);
                }
                else
                {
                    ++index;
                    continue;
                }

                LOGGER.debug("adding assembly linked pair({}) to chain({})", assemblyLink.toString(), currentChain.getId());
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

            SvLinkedPair closestStartPair = findNextLinkedPair(requiredLinks, unlinkedSvList, chainFirstSV, chainFirstUnlinkedOnStart);
            SvLinkedPair closestLastPair = findNextLinkedPair(requiredLinks, unlinkedSvList, chainLastSV, chainLastUnlinkedOnStart);

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
                        LOGGER.debug("adding linked pair({}) to chain({}) start with length({})",
                                closestStartPair.toString(), currentChain.getId(), closestStartPair.length());
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
                        LOGGER.debug("adding linked pair({}) to chain({}) end with length({})",
                                closestLastPair.toString(), currentChain.getId(), closestLastPair.length());
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
                    LOGGER.debug("sample({}) cluster({}) adding {} chain({}) with {} linked pairs:",
                            mSampleId, mCluster.getId(), chainComplete ? "complete" : "partial", currentChain.getId(), currentChain.getLinkCount());
                    currentChain.logLinks();
                }

                if(chainComplete)
                    break;
            }
        }
    }

    private SvLinkedPair findNextLinkedPair(final List<SvLinkedPair> chainedPairs, final List<SvVarData> svList, final SvVarData var, boolean useStart)
    {
        if(var.isAssemblyMatched(useStart)) // this breakend will be added when the assembly link is added
            return null;

        // find the shortest templated insertion which links to this variant
        SvLinkedPair newPair = null;

        for(final SvVarData otherVar : svList)
        {
            if(otherVar.type() == INS)
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

                if (tiLength <= 0 || (newPair != null && tiLength >= newPair.length()))
                    continue;

                // check if already in a linked pair
                boolean isValid = true;
                for(final SvLinkedPair pair : chainedPairs)
                {
                    if(pair.hasVariantBE(otherVar, otherVarStart))
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
                    // the breakends used in the link match each other
                    for(final SvLinkedPair pair : chainedPairs)
                    {
                        // check that if this link has the same SVs as another link due to replicated SVs, that
                        // the breakends used in the link match each other
                        if(testPair.sameVariants(pair) && !testPair.matches(pair, true))
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

    private static boolean linkChains(List<SvChain> chains, int reqSvCount)
    {
        // reconcile chains together if possible
        if (chains.size() <= 1)
            return false;

        int chainsRemoved = 0;

        for (int firstIndex = 0; firstIndex < chains.size(); ++firstIndex)
        {
            final SvChain firstChain = chains.get(firstIndex);

            if (firstChain.getSvCount() == reqSvCount || firstChain.isClosedLoop())
                continue;

            int nextIndex = firstIndex + 1;
            while (nextIndex < chains.size())
            {
                final SvChain nextChain = chains.get(nextIndex);

                if (nextChain.getSvCount() == reqSvCount || nextChain.isClosedLoop())
                {
                    ++nextIndex;
                    continue;
                }

                if (reconcileChains(firstChain, nextChain))
                {
                    chains.remove(nextIndex);
                    ++chainsRemoved;
                }
                else
                {
                    ++nextIndex;
                }
            }
        }

        return chainsRemoved > 0;
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

    private boolean findCompletedLinks(
            final List<SvLinkedPair> requiredLinks, final List<SvLinkedPair> existingLinks,
            final List<SvChain> partialChains, List<SvLinkedPair> testLinks, int currentIndex)
    {
        if (currentIndex >= testLinks.size())
            return false;

        if(mCompleteChains.size() >= MAX_COMPLETE_CHAINS)
            return false;

        if(mPathInterations >= MAX_PATH_ITERATIONS)
            return false;

        ++mPathInterations;

        if(mLogPathFinding)
        {
            LOGGER.debug("currentIndex({}) links(existing={} test={}) chains(full={} part={}) partials({}) reqSvCount({}) pathIters({})",
                    currentIndex, existingLinks.size(), testLinks.size(), mCompleteChains.size(), mIncompleteChains.size(), partialChains.size(), mReqChainCount, mPathInterations);
        }

        boolean linksAdded = false;

        // try adding each of the remaining links recursively to get a set of completely-links breakends
        // and if so, add this collection of potential chains
        for (int i = currentIndex; i < testLinks.size(); ++i)
        {
            final SvLinkedPair testLink = testLinks.get(i);

            // check this ne test link against both the existing links and the existing chains
            if (!canAddLinkedPair(existingLinks, partialChains, testLink, mReqChainCount))
                continue;

            linksAdded = true;

            // must be a deep copy to avoid changing the set of partial chains at this level of the search
            List<SvChain> workingChains = Lists.newArrayList();
            for (final SvChain chain : partialChains)
            {
                workingChains.add(new SvChain(chain));
            }

            // add this link to any chains if possible
            addLinkToChains(testLink, workingChains);

            // every new SV is now part of a partial chain (even if only a single link)
            List<SvLinkedPair> workingLinks = Lists.newArrayList();
            workingLinks.addAll(existingLinks);
            workingLinks.add(testLink);

            // reconcile chains together if possible
            linkChains(workingChains, mReqChainCount);

            // and check whether any chains are complete
            int chainIndex = 0;
            while (chainIndex < workingChains.size())
            {
                final SvChain chain = workingChains.get(chainIndex);

                if (chain.getSvCount() < mReqChainCount)
                {
                    ++chainIndex;
                    continue;
                }

                boolean skipChain = false;
                if(mCluster.hasReplicatedSVs() && !validateReplicatingChain(chain))
                {
                    skipChain = true;
                }
                else if(mHasExistingChains && !validatesCopyNumberData(chain))
                {
                    skipChain = true;
                }
                else if(isDuplicateChain(chain, mCompleteChains))
                {
                    // LOGGER.debug("skipping duplicate complete chain:");
                    // chain.logLinks();

                    skipChain = true;
                }

                boolean hasRequiredLinks = chain.hasLinks(requiredLinks);

                if (!skipChain && hasRequiredLinks)
                {
                    chain.setId(mCompleteChains.size());
                    mCompleteChains.add(chain);

                    if(mLogVerbose)
                    {
                        LOGGER.debug("iters({} of {} current={}) existingLinks({}) completedChains({}) workingChains({}) - found complete chain({} svs={}):",
                                i, testLinks.size(), currentIndex, existingLinks.size(), mCompleteChains.size(), workingChains.size(),
                                chain.getId(), chain.getSvCount());

                        chain.logLinks();
                    }
                }

                workingChains.remove(chainIndex);
            }

            if (workingChains.isEmpty())
            {
                // this search path is done so no need to keep searching for more links to add
                break;
            }

            if(mCompleteChains.size() >= MAX_COMPLETE_CHAINS)
                return false;

            if(mLogPathFinding)
            {
                LOGGER.debug("iters({} of {} current={}) existingLinks({}) chains(full={} part={}) workingChains({}) - continuing search",
                        i, testLinks.size(), currentIndex, existingLinks.size(), mCompleteChains.size(), mIncompleteChains.size(), workingChains.size());
            }

            // continue the search, moving on to try adding the next test link
            boolean hasMoreTestLinks = findCompletedLinks(requiredLinks, workingLinks, workingChains, testLinks, i + 1);
            boolean hasCompleteChains = hasCompleteChains(mReqChainCount);

            if (!hasMoreTestLinks && !hasCompleteChains && !mRequireFullChains)
            {
                // add any sufficiently long chains since this search path has been exhausted
                for (final SvChain chain : workingChains)
                {
                    if (chain.getSvCount() < mMinIncompleteChainCount)
                        continue;

                    // check for duplicates
                    if(isDuplicateChain(chain, mIncompleteChains) || shorterThanExistingIncompletes(chain))
                        continue;

                    chain.setId(mIncompleteChains.size());
                    mIncompleteChains.add(chain);

                    if(mLogVerbose)
                    {
                        LOGGER.debug("iters({} of {} current={}) existingLinks({}) chains(full={} part={}) workingChains({}) - found incomplete chain({} svs={})",
                                i, testLinks.size(), currentIndex, existingLinks.size(), mCompleteChains.size(), mIncompleteChains.size(), workingChains.size(),
                                chain.getId(), chain.getSvCount());
                    }
                }
            }
        }

        return linksAdded;
    }

    private boolean isDuplicateChain(final SvChain chain, final List<SvChain> chains)
    {
        for (final SvChain existingChain : chains)
        {
            if (chain.isIdentical(existingChain))
                return true;
        }

        return false;
    }

    private boolean shorterThanExistingIncompletes(final SvChain chain)
    {
        for (final SvChain existingChain : mIncompleteChains)
        {
            if (chain.getSvCount() >= existingChain.getSvCount())
                continue;

            // if this new chain is shorter and has some of the same links, then skip if
            for(final SvVarData var : chain.getSvList())
            {
                if(existingChain.getSvList().contains(var))
                {
                    return true;
                }
            }
        }

        return false;
    }

    private boolean hasCompleteChains(int reqSvCount)
    {
        for (final SvChain chain : mCompleteChains)
        {
            if (chain.getSvCount() == reqSvCount)
                return true;
        }

        return false;
    }

    private static void addLinkToChains(SvLinkedPair linkedPair, List<SvChain> chains)
    {
        // add to an existing chain or create a new one
        for(final SvChain chain : chains)
        {
            if(chain.hasLinkedPair(linkedPair))
                return;

            if(chain.hasLinkClash(linkedPair))
                continue;

            if(chain.canAddLinkedPairToStart(linkedPair))
            {
                chain.addLink(linkedPair, true);
                return;
            }

            if(chain.canAddLinkedPairToEnd(linkedPair))
            {
                chain.addLink(linkedPair, false);
                return;
            }
        }

        SvChain newChain = new SvChain(0);
        newChain.addLink(linkedPair, true);
        chains.add(newChain);
    }

    private static boolean canAddLinkedPair(final List<SvLinkedPair> existingLinks, final List<SvChain> partialChains, final SvLinkedPair testLink, int requiredSvCount)
    {
        // first check that the breakends aren't already used
        for(SvLinkedPair linkedPair : existingLinks)
        {
            if (linkedPair.hasLinkClash(testLink))
                return false;

            // check that if this link has the same SVs as another link due to replicated SVs, that
            // the breakends used in the link match each other
            if(testLink.sameVariants(linkedPair) && !testLink.matches(linkedPair, true))
                return false;
        }

        // then check that this linked pair doesn't close a chain of a smaller size than one involving all required SVs
        for(final SvChain chain : partialChains)
        {
            boolean linkCouldCloseChain = chain.linkWouldCloseChain(testLink);

            if(!linkCouldCloseChain)
                continue;

            // work out whether all SVs would be accounted for
            int uniqueSVCount = chain.getSvCount();

            if(!chain.getSvList().contains(testLink.first()))
                ++uniqueSVCount;

            if(!chain.getSvList().contains(testLink.second()))
                ++uniqueSVCount;

            if(uniqueSVCount < requiredSvCount)
            {
                // closing this chain would make it too short
                return false;
            }
            else
            {
                // this would be a looped chain involving all required SVs
                return true;
            }
        }

        return true;
    }

    private boolean validateReplicatingChain(final SvChain chain)
    {
        // any chain which exhibits copy number growth (eg BFB) must start or end on the outermost INV facing back towards the centromere
        for(int i = 0; i < 2; ++i)
        {
            final SvVarData unlinkedSv = (i == 0) ? chain.getFirstSV() : chain.getLastSV();

            if(unlinkedSv.type() != INV && !unlinkedSv.isTranslocation())
                continue;

            boolean useStartPosition = unlinkedSv.orientation(true) == -1;
            // String arm = unlinkedSv.orientation(true) == 1 ? unlinkedSv.arm(false) : unlinkedSv.arm(true);
            long chainOuterPosition = unlinkedSv.position(useStartPosition);

            // check for another INV further out towards the telomere than this one
            boolean hasFurtherOutSV = false;
            for(SvVarData var : chain.getSvList())
            {
                if(var.isNullBreakend() && !useStartPosition)
                    continue;

                if((useStartPosition && var.position(useStartPosition) < chainOuterPosition)
                        || (!useStartPosition && var.position(useStartPosition) > chainOuterPosition))
                {
                    if(mLogVerbose)
                    {
                        LOGGER.debug("chain({}) {} SV({} pos={}) not outermost vs SV({} pos={})",
                                chain.getId(), (i==0)  ? "start" : "end", unlinkedSv.id(), chainOuterPosition, var.id(), var.position(useStartPosition));
                    }

                    hasFurtherOutSV = true;
                    break;
                }
            }

            if(!hasFurtherOutSV)
                return true;
        }

        return false;
    }


    private boolean validatesCopyNumberData(final SvChain chain)
    {
        /*
        final Map<String, List<SvCNData>> chrCNDataMap = mCluster.getChrCNData();
        Map<String, List<Integer>> chrSegmentPathMap = new HashMap();

        final List<SvLinkedPair> linkedPairs = chain.getLinkedPairs();
        final List<SvVarData> svList = chain.getSvList();

        // walk through the chain looking at each segment formed either by an SV or the templated insertion from 2 SVs
        // find the corresponding copy-number segment and register a path through (literally adding +1 to the path)
        for(int i = 0; i < linkedPairs.size(); ++i)
        {
            final SvLinkedPair pair = linkedPairs.get(i);

            String chromosome = "";
            long startPosition = 0;
            long endPosition = 0;

            if(i == 0 || i == linkedPairs.size()-1)
            {
                // the first open link on the chain is taken as if coming from either the telomere or centromere
                // for the last open SV also register its path as heading out to telomere/centromere
                final SvVarData var = (i == 0) ? svList.get(i) : svList.get(i+1);;

                boolean useStart = (var == pair.first() && pair.firstUnlinkedOnStart()) || (var == pair.second() && pair.secondUnlinkedOnStart());

                if(!var.isNullBreakend() || useStart)
                {
                    chromosome = var.chromosome(useStart);

                    byte orientation = var.orientation(useStart);

                    if (orientation == 1)
                    {
                        startPosition = 0;
                        endPosition = var.position(useStart);
                    }
                    else
                    {
                        startPosition = var.position(useStart);
                        endPosition = mUtils.getChromosomalArmLength(chromosome, var.arm(useStart));
                    }

                    registerCopyNumberSegment(chrCNDataMap, chrSegmentPathMap, chromosome, startPosition, endPosition);
                }
            }

            // subsequent SVs path a back between the linked ends along the templated insertion length
            long firstPos = pair.first().position(pair.firstLinkOnStart());
            long secondPos = pair.second().position(pair.secondLinkOnStart());
            startPosition = firstPos < secondPos ? firstPos : secondPos;
            endPosition = secondPos > firstPos ? secondPos : firstPos;
            chromosome = pair.first().chromosome(pair.firstLinkOnStart());

            registerCopyNumberSegment(chrCNDataMap, chrSegmentPathMap, chromosome, startPosition, endPosition);
        }

        // compare the path segment counts with the actual copy number profile
        boolean cnProfileMatched = true;
        for(Map.Entry<String, List<Integer>> entry : chrSegmentPathMap.entrySet())
        {
            final String chromosome = entry.getKey();
            final List<SvCNData> copyNumberData = chrCNDataMap.get(chromosome);
            List<Integer> segmentPathList = chrSegmentPathMap.get(chromosome);

            for(int i = 0; i < segmentPathList.size(); ++i)
            {
                final SvCNData cnData = copyNumberData.get(i);
                int pathSegmentCount = segmentPathList.get(i);

                LOGGER.debug(String.format("cluster(%d) chromosome(%s) seg %d: position(%d) copyNumber(%.1f) pathCount(%d)",
                        mCluster.getId(), chromosome, i, cnData.startPos(), cnData.copyNumber(), pathSegmentCount));

                if(cnData.copyNumber() != pathSegmentCount)
                {
                    cnProfileMatched = false;
                }
            }
        }
        */

        // return cnProfileMatched;
        return true;
    }

    private void registerCopyNumberSegment(
            final Map<String, List<SvCNData>> chrCNDataMap, Map<String, List<Integer>> chrSegmentPathMap,
            final String chromosome, long startPosition, long endPosition)
    {
        if(!chrCNDataMap.containsKey(chromosome))
        {
            LOGGER.error("sample({}) cluster({}) chromosome({}) not found with positions({} -> {})",
                    mSampleId, mCluster.getId(), chromosome, startPosition, endPosition);
            return;
        }

        if(startPosition >= endPosition)
        {
            LOGGER.error("sample({}) cluster({}) invalid segment chr({}) positions({} -> {})",
                    mSampleId, mCluster.getId(), chromosome, startPosition, endPosition);
            return;
        }

        final List<SvCNData> cnList = chrCNDataMap.get(chromosome);

        List<Integer> segmentPathList = chrSegmentPathMap.get(chromosome);
        if(segmentPathList == null)
        {
            segmentPathList = Lists.newArrayList();

            for(int i = 0; i < cnList.size(); ++i)
                segmentPathList.add(i, 0);

            chrSegmentPathMap.put(chromosome, segmentPathList);
        }

        // coming in from the right, so need to mark off all segments going right from this point
        boolean foundStart = false;

        for(int i = 0; i < cnList.size(); ++i)
        {
            final SvCNData cnData = cnList.get(i);

            if(!foundStart)
            {
                if(cnData.startPos() != startPosition)
                    continue;

                foundStart = true;
            }
            else
            {
                // check if have proceeded too far
                if(cnData.startPos() > endPosition)
                    break;
            }

            int segCount = segmentPathList.get(i);
            segmentPathList.set(i, segCount+ 1);
        }

        if(!foundStart)
        {
            LOGGER.error("sample({}) cluster({}) segment not found for chr({}) positions({} -> {})",
                    mSampleId, mCluster.getId(), chromosome, startPosition, endPosition);
        }
    }

}
