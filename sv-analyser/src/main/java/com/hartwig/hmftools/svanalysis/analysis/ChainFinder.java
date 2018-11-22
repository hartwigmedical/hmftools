package com.hartwig.hmftools.svanalysis.analysis;

import static java.lang.Math.abs;
import static java.lang.Math.max;
import static java.lang.Math.min;
import static java.lang.Math.round;
import static java.lang.Math.sin;

import static com.hartwig.hmftools.common.variant.structural.StructuralVariantType.INS;
import static com.hartwig.hmftools.common.variant.structural.StructuralVariantType.INV;
import static com.hartwig.hmftools.common.variant.structural.StructuralVariantType.SGL;
import static com.hartwig.hmftools.svanalysis.analysis.LinkFinder.MIN_TEMPLATED_INSERTION_LENGTH;
import static com.hartwig.hmftools.svanalysis.analysis.LinkFinder.areLinkedSection;
import static com.hartwig.hmftools.svanalysis.analysis.SvUtilities.calcTypeCount;
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
    private static double LOW_PLOIDY_LIMIT = 0.05;

    private static final Logger LOGGER = LogManager.getLogger(ChainFinder.class);

    final SvUtilities mUtils;
    private String mSampleId;
    private SvCluster mCluster;
    private List<SvLinkedPair> mAssemblyLinkedPairs;
    private List<SvLinkedPair> mInferredLinkedPairs;
    private List<SvChain> mCompleteChains;
    private List<SvChain> mIncompleteChains;
    private boolean mHasExistingChains;
    private int mMinIncompleteChainCount;
    private int mReqChainCount;
    private boolean mLogVerbose;

    PerformanceCounter mContinuousFinderPc;

    public ChainFinder(final SvUtilities utils)
    {
        mUtils = utils;
        mCompleteChains = Lists.newArrayList();
        mIncompleteChains = Lists.newArrayList();
        mLogVerbose = false;
        mReqChainCount = 0;
        mMinIncompleteChainCount = 0;

        mContinuousFinderPc = new PerformanceCounter("Continuous");
    }

    public final PerformanceCounter getContinuousFinderPc() { return mContinuousFinderPc; }

    public void initialise(final String sampleId, SvCluster cluster)
    {
        mSampleId = sampleId;
        mCluster = cluster;
        mHasExistingChains = !mCluster.getChains().isEmpty();
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

    public void setLogVerbose(boolean toggle) { mLogVerbose = toggle; }

    public final List<SvChain> getCompleteChains() { return mCompleteChains; }

    public boolean formClusterChains()
    {
        // take the assembly links as a given and then try out the inferred links to see if a single chain can be formed from all the breakends

        // only factor in templated insertions to form chains, even if they also exhibit DBs
        List<SvLinkedPair> inferredTIs = collectRelevantLinkedPairs();

        List<SvVarData> svList = collectRelevantSVs();

        mReqChainCount = svList.size();

        /*
        // check whether there are enough potential links to form full or near-to-full chains
        if (!mCluster.hasSubClusters() && !checkLinksPotential(svList, inferredTIs))
            return false;
        */

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
        return mInferredLinkedPairs;

        /*
        List<SvLinkedPair> inferredTIs = Lists.newArrayList();
        for (SvLinkedPair pair : mInferredLinkedPairs)
        {
            if (!(pair.linkType() == LINK_TYPE_TI || pair.linkType() == LINK_TYPE_SGL))
                continue;

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

            inferredTIs.add(pair);
        }

        return inferredTIs;
        */
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

            int maxSubChainCount = 0;
            for(final SvChain chain : mCluster.getChains())
            {
                LOGGER.debug("cluster({}) has existing chain({} svs={})", mCluster.getId(), chain.getId(), chain.getSvCount());

                maxSubChainCount = max(maxSubChainCount, chain.getSvCount());
            }

            mMinIncompleteChainCount = max(maxSubChainCount, mMinIncompleteChainCount);
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

        // maintain links and SVs already allocated to partial chains, but add the open ends of chains to the unlinked list
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
            if(!chainLinkAdded || currentChain == null)
            {
                if(!partialChains.isEmpty())
                {
                    currentChain = new SvChain(partialChains.get(0));
                    partialChains.remove(0);

                    if (mLogVerbose)
                    {
                        LOGGER.debug("sample({}) cluster({}) building from existing chain({}) with {} SVs",
                                mSampleId, mCluster.getId(), currentChain.getId(), currentChain.getSvCount());
                    }

                    for(final SvLinkedPair pair : currentChain.getLinkedPairs())
                    {
                        if(!chainedPairs.contains(pair))
                            chainedPairs.add(pair);

                        reduceRemainingLists(pair, unlinkedSvList, remainingStartLinks);
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
                        assemblyLink.toString(), currentChain.getId(), addToStart ? "start" : "end");

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
                        LOGGER.debug("adding linked pair({}) to chain({}) start",
                                closestStartPair.toString(), currentChain.getId());
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
                                closestLastPair.toString(), currentChain.getId());
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

                currentChain = null;

                if(chainComplete)
                    break;
            }
        }

        if(!chainComplete && !chainsList.isEmpty() && unlinkedSvList.size() == 1 && unlinkedSvList.get(0).type() == SGL)
        {
            tryUnlinkedSingle(chainsList, unlinkedSvList.get(0));
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

                if (tiLength <= MIN_TEMPLATED_INSERTION_LENGTH || (newPair != null && tiLength >= newPair.length()))
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
                LOGGER.debug("adding linked pair({}) with SGL to chain({}) start", newPair.toString(), chain.getId());
                amendedChain = chain;
                addToStart = true;
                break;
            }

            if (chain.getLastSV().type() == SGL && !chain.getLastSV().equals(singleSV, true))
            {
                newPair = new SvLinkedPair(chain.getLastSV(), singleSV, LINK_TYPE_TI, false, false);
                LOGGER.debug("adding linked pair({}) with SGL to chain({}) end", newPair.toString(), chain.getId());
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
                    LOGGER.debug("adding existing chain({}) to current chain({})", chain.getId(), amendedChain.getId());
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
                mCluster.getId(), chain.getId(),
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
