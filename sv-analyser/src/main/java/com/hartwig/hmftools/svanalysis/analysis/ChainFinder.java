package com.hartwig.hmftools.svanalysis.analysis;

import static java.lang.Math.abs;
import static java.lang.Math.max;
import static java.lang.Math.min;
import static java.lang.Math.round;

import static com.hartwig.hmftools.common.variant.structural.StructuralVariantType.INS;
import static com.hartwig.hmftools.svanalysis.analysis.ClusterAnalyser.MAX_COPY_NUMBER_DIFF;
import static com.hartwig.hmftools.svanalysis.analysis.ClusterAnalyser.MAX_COPY_NUMBER_DIFF_PERC;
import static com.hartwig.hmftools.svanalysis.types.SvClusterData.RELATION_TYPE_NEIGHBOUR;
import static com.hartwig.hmftools.svanalysis.types.SvLinkedPair.LINK_TYPE_SGL;
import static com.hartwig.hmftools.svanalysis.types.SvLinkedPair.LINK_TYPE_TI;

import java.util.List;
import java.util.Map;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.svanalysis.types.SvCNData;
import com.hartwig.hmftools.svanalysis.types.SvChain;
import com.hartwig.hmftools.svanalysis.types.SvClusterData;
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
    private int mReqSvCount;
    private int mPathInterations;
    private boolean mLogVerbose;

    public ChainFinder(final SvUtilities utils)
    {
        mUtils = utils;
        mRequireFullChains = false;
        mCompleteChains = Lists.newArrayList();
        mIncompleteChains = Lists.newArrayList();
        mLogVerbose = false;
        mPathInterations = 0;
        mReqSvCount = 0;
        mMinIncompleteChainCount = 0;
    }

    public void initialise(final String sampleId, SvCluster cluster)
    {
        mSampleId = sampleId;
        mCluster = cluster;
        mHasExistingChains = !mCluster.getChains().isEmpty();
        mPathInterations = 0;
        mReqSvCount = 0;
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

        List<SvClusterData> svList = collectRelevantSVs();

        mReqSvCount = svList.size();

        if(!mRequireFullChains)
            mMinIncompleteChainCount = (int)round(MIN_CHAIN_PERCENT * mReqSvCount);

        // check whether there are enough potential links to form full or near-to-full chains
        if (!checkLinksPotential(svList, inferredTIs))
            return false;

        // for now if there are too many potential inferred links, cull the list first
        /*
        int reqInferredLinks = svList.size() - mAssemblyLinkedPairs.size();
        if (reqInferredLinks > 10 && inferredTIs.size() > 10)
        {
            reduceInferredToShortestLinks(inferredTIs, mCompleteChains);

            if (!checkLinksPotential(svList, inferredTIs))
                return false;
        }
        */

        if (mCluster.getCount() >= 4)
        {
            LOGGER.debug("sample({}) cluster({}) assemblyLinks({}) inferredLinks({} of total={}) svCount({} all={})) existingChains({})",
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

    private List<SvClusterData> collectRelevantSVs()
    {
        List<SvClusterData> svList = Lists.newArrayList();

            // remove any non-applicable or suspect SVs
        for(final SvClusterData var : mCluster.getSVs())
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

    private boolean checkLinksPotential(List<SvClusterData> svList, List<SvLinkedPair> inferredLinkedPairs)
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

        for (final SvClusterData var : svList)
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

    private boolean formChains(List<SvClusterData> svList, List<SvLinkedPair> inferredLinkedPairs)
    {
        // find all the combinations of linked pairs where every breakend end except at most 2 are covered by a linked pair
        int clusterSvCount = svList.size();

        List<SvLinkedPair> trialLinkedPairs = Lists.newArrayList();
        trialLinkedPairs.addAll(mAssemblyLinkedPairs);

        List<SvChain> partialChains = Lists.newArrayList();

        if(mHasExistingChains)
        {
            partialChains.addAll(mCluster.getChains());
            linkChains(partialChains, clusterSvCount);

            if(!mRequireFullChains)
            {
                int maxSubChainCount = 0;
                for(final SvChain chain : mCluster.getChains())
                {
                    maxSubChainCount = max(maxSubChainCount, chain.getSvCount());
                }

                mMinIncompleteChainCount = max(maxSubChainCount, mMinIncompleteChainCount);
            }
        }

        // first form any partial chains out of the assembly linked pairs
        for (SvLinkedPair linkedPair : mAssemblyLinkedPairs)
        {
            addLinkToChains(linkedPair, partialChains);
        }

        if(partialChains.size() == 1 && partialChains.get(0).getSvCount() == clusterSvCount)
        {
            // the assemblies formed a single chain with all the required link
            SvChain completeChain = partialChains.get(0);
            completeChain.setId(mCompleteChains.size());
            mCompleteChains.add(completeChain);
            partialChains.clear();

            LOGGER.debug("found complete chain({} svs={}) from assembly links",
                    completeChain.getId(), completeChain.getSvCount());

            completeChain.logLinks();
        }
        else
        {
            // check if there are too many links to check
            /*
            int reqInferredLinks = svList.size() - mAssemblyLinkedPairs.size() - 1;
            int testLinkCount = inferredLinkedPairs.size();

            if(reqInferredLinks > 20 || (reqInferredLinks >= 10 && testLinkCount >= 40))
                return false;

            if(testLinkCount > 10 && testLinkCount > reqInferredLinks)
            {
                long calcCount = combination(testLinkCount, reqInferredLinks);

                if (calcCount > 10000)
                    return false;
            }
            */

            findCompletedLinks(clusterSvCount, mAssemblyLinkedPairs, mAssemblyLinkedPairs, partialChains, inferredLinkedPairs, 0);

            if(mPathInterations >= MAX_PATH_ITERATIONS)
            {
                LOGGER.info("sample({}) cluster({}) max iterations reached: chains(full={} part={}) svCount({}) links(assembly={} infer={})",
                        mSampleId, mCluster.getId(), mCompleteChains.size(), mIncompleteChains.size(),
                        svList.size(), mAssemblyLinkedPairs.size(), inferredLinkedPairs.size());
            }
        }

        // first check for any complete chains, and if found, add the shortest one
        if(cacheCompleteChain(clusterSvCount))
            return true;

        // otherwise add the longest mutually exclusive chains
        cacheIncompleteChains();

        return false;
    }

    private boolean cacheCompleteChain(int reqChainCount)
    {
        if(mCompleteChains.isEmpty())
            return false;

        // cache any complete chain found
        int shortestLength = -1;
        SvChain shortestFullChain = null;
        for (SvChain chain : mCompleteChains)
        {
            chain.recalcLength();

            if (chain.getLinkCount() < 2 || chain.getSvCount() != reqChainCount)
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

        List<SvClusterData> chainedSVs = Lists.newArrayList();

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
                    for (final SvClusterData var : chain.getSvList())
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
                LOGGER.info("sample({}) cluster({}) adding incomplete chain({}) length({}) with {} linked pairs",
                        mSampleId, mCluster.getId(), maxLengthChain.getId(), maxLengthChain.getLength(), maxLengthChain.getLinkCount());

                maxLengthChain.logLinks();
            }

            mCluster.addChain(maxLengthChain);
            chainedSVs.addAll(maxLengthChain.getSvList());
            mIncompleteChains.remove(maxLengthChain);
        }
    }

    private boolean findCompletedLinks(
            int reqSvCount, final List<SvLinkedPair> requiredLinks, final List<SvLinkedPair> existingLinks,
            final List<SvChain> partialChains, List<SvLinkedPair> testLinks, int currentIndex)
    {
        if (currentIndex >= testLinks.size())
            return false;

        if(mCompleteChains.size() >= MAX_COMPLETE_CHAINS)
            return false;

        if(mPathInterations >= MAX_PATH_ITERATIONS)
            return false;

        ++mPathInterations;

        if(mLogVerbose)
        {
            LOGGER.debug("currentIndex({}) links(existing={} test={}) chains(full={} part={}) partials({}) reqSvCount({}) pathIters({})",
                    currentIndex, existingLinks.size(), testLinks.size(), mCompleteChains.size(), mIncompleteChains.size(), partialChains.size(), reqSvCount, mPathInterations);
        }

        boolean linksAdded = false;

        // try adding each of the remaining links recursively to get a set of completely-links breakends
        // and if so, add this collection of potential chains
        for (int i = currentIndex; i < testLinks.size(); ++i)
        {
            final SvLinkedPair testLink = testLinks.get(i);

            // check this ne test link against both the existing links and the existing chains
            if (!canAddLinkedPair(existingLinks, partialChains, testLink, reqSvCount))
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
            linkChains(workingChains, reqSvCount);

            // and check whether any chains are complete
            int chainIndex = 0;
            while (chainIndex < workingChains.size())
            {
                final SvChain chain = workingChains.get(chainIndex);

                if (chain.getSvCount() < reqSvCount)
                {
                    ++chainIndex;
                    continue;
                }

                // check if validates the copy number profile
                if(mHasExistingChains && !validatesCopyNumberData(chain))
                    continue;

                boolean isDuplicate = isDuplicateChain(chain, mCompleteChains);

                /*
                if(isDuplicate)
                {
                    LOGGER.debug("skipping duplicate complete chain:");
                    chain.logLinks();
                }
                */

                boolean hasRequiredLinks = chain.hasLinks(requiredLinks);

                if (!isDuplicate && hasRequiredLinks)
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

            if(mLogVerbose)
            {
                LOGGER.debug("iters({} of {} current={}) existingLinks({}) chains(full={} part={}) workingChains({}) - continuing search",
                        i, testLinks.size(), currentIndex, existingLinks.size(), mCompleteChains.size(), mIncompleteChains.size(), workingChains.size());
            }

            // continue the search, moving on to try adding the next test link
            boolean hasMoreTestLinks = findCompletedLinks(reqSvCount, requiredLinks, workingLinks, workingChains, testLinks, i + 1);
            boolean hasCompleteChains = hasCompleteChains(reqSvCount);

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
            for(final SvClusterData var : chain.getSvList())
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

    private boolean validatesCopyNumberData(final SvChain chain)
    {
        /*
        final Map<String, List<SvCNData>> chrCNDataMap = mCluster.getChrCNData();
        Map<String, List<Integer>> chrSegmentPathMap = new HashMap();

        final List<SvLinkedPair> linkedPairs = chain.getLinkedPairs();
        final List<SvClusterData> svList = chain.getSvList();

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
                final SvClusterData var = (i == 0) ? svList.get(i) : svList.get(i+1);;

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

    // old methods which don't use assembly info
    public void findSvChains(final String sampleId, SvCluster cluster)
    {
        if(!cluster.getChains().isEmpty())
            return;

        // findSvChains(sampleId, cluster, cluster.getLinkedPairs(), cluster.getChains());
    }

    public void findSvChains(final String sampleId, SvCluster cluster, List<SvLinkedPair> linkedPairs, List<SvChain> chainsList)
    {
        if(linkedPairs.isEmpty() || linkedPairs.size() < 2)
            return;

        List<SvLinkedPair> workingLinkedPairs = Lists.newArrayList();
        workingLinkedPairs.addAll(linkedPairs);

        LOGGER.debug("cluster({}) attempting to find chained SVs from {} linked pairs", cluster.getId(), workingLinkedPairs.size());

        while(workingLinkedPairs.size() >= 2)
        {
            // start with a single linked pair
            // for each of its ends (where the first BE is labelled 'first', and the second labelled 'last'),
            // search for the closest possible linking BE from another linked pair
            // for BEs to link they must be facing (like a TI)

            SvLinkedPair linkedPair = workingLinkedPairs.get(0);
            workingLinkedPairs.remove(0);

            int chainId = chainsList.size()+1;
            SvChain currentChain = new SvChain(chainId);

            LOGGER.debug("sample({}) cluster({}) starting chain({}) with linked pair({})",
                    sampleId, cluster.getId(), chainId, linkedPair.toString());

            currentChain.addLink(linkedPair, true);

            while(!workingLinkedPairs.isEmpty())
            {
                // now search the remaining SVs for links at either end of the current chain
                SvClusterData beFirst = currentChain.getFirstSV();
                boolean chainFirstUnlinkedOnStart = currentChain.firstLinkOpenOnStart();
                SvClusterData beLast = currentChain.getLastSV();
                boolean chainLastUnlinkedOnStart = currentChain.lastLinkOpenOnStart();

                SvLinkedPair closestStartPair = null;
                int closestStartLen = -1;

                SvLinkedPair closestLastPair = null;
                int closestLastLen = -1;

                for(SvLinkedPair pair : workingLinkedPairs) {

                    // first check for a linked pair which has the same variant (but the other BE) to the unlinked on
                    if((beFirst.equals(pair.first()) && chainFirstUnlinkedOnStart == pair.firstLinkOnStart())
                            || (beFirst.equals(pair.second()) && chainFirstUnlinkedOnStart == pair.secondLinkOnStart()))
                    {
                        closestStartPair = pair;
                        closestStartLen = 0; // to prevent another match
                    }

                    if((beLast.equals(pair.first()) && chainLastUnlinkedOnStart == pair.firstLinkOnStart())
                            || (beLast.equals(pair.second()) && chainLastUnlinkedOnStart == pair.secondLinkOnStart()))
                    {
                        closestLastPair = pair;
                        closestLastLen = 0; // to prevent another match
                    }
                }

                if(closestStartPair == null && closestLastPair == null)
                {
                    break;
                }

                if(closestStartPair != null)
                {
                    LOGGER.debug("adding linked pair({}) on chain start({}) with length({})",
                            closestStartPair.toString(), beFirst.posId(chainFirstUnlinkedOnStart), closestStartLen);

                    // add this to the chain at the start
                    currentChain.addLink(closestStartPair, true);
                    workingLinkedPairs.remove(closestStartPair);
                }

                if(closestLastPair != null & closestLastPair != closestStartPair)
                {
                    LOGGER.debug("adding linked pair({}) on chain end({}) with length({})",
                            closestLastPair.toString(), beLast.posId(chainLastUnlinkedOnStart), closestLastLen);

                    // add this to the chain at the start
                    currentChain.addLink(closestLastPair, false);
                    workingLinkedPairs.remove(closestLastPair);
                }
            }

            if(currentChain.getLinkCount() > 1)
            {
                LOGGER.debug("sample({}) cluster({}) adding chain({}) with {} linked pairs:",
                        sampleId, cluster.getId(), currentChain.getId(), currentChain.getLinkCount());

                chainsList.add(currentChain);

                for(int i = 0; i < currentChain.getLinkCount(); ++i)
                {
                    final SvLinkedPair pair = currentChain.getLinkedPairs().get(i);
                    LOGGER.debug("sample({}) cluster({}) chain({}) {}: pair({}) {} len={}",
                            sampleId, cluster.getId(), currentChain.getId(), i, pair.toString(), pair.linkType(), pair.length());
                }
            }
        }
    }


}
