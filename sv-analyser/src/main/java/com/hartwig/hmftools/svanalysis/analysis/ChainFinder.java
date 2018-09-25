package com.hartwig.hmftools.svanalysis.analysis;

import static com.hartwig.hmftools.svanalysis.types.SvLinkedPair.hasCompleteLinkedPairsList;

import java.util.List;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.svanalysis.types.SvChain;
import com.hartwig.hmftools.svanalysis.types.SvClusterData;
import com.hartwig.hmftools.svanalysis.types.SvLinkedPair;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;

public class ChainFinder
{
    private static final Logger LOGGER = LogManager.getLogger(ChainFinder.class);

    public static boolean assessClusterChaining(final String sampleId, SvCluster cluster, List<SvLinkedPair> assemblyLinkedPairs, List<SvLinkedPair> inferredLinkedPairs)
    {
        // take the assembly links as a given and then try out the inferred links to see if a single chain can be formed from all the breakends

        // find all the combinations of linked pairs where every breakend end except at most 2 are covered by a linked pair
        List<SvLinkedPair> trialLinkedPairs = Lists.newArrayList();
        trialLinkedPairs.addAll(assemblyLinkedPairs);

        List<SvChain> completeChains = Lists.newArrayList();
        List<SvChain> partialChains = Lists.newArrayList();

        // first form any partial chains out of the assembly linked pairs
        for(SvLinkedPair linkedPair : assemblyLinkedPairs)
        {
            addLinkToChains(linkedPair, partialChains);
        }

        findCompletedLinks(cluster.getSVs(), assemblyLinkedPairs, assemblyLinkedPairs, partialChains, inferredLinkedPairs, 0, completeChains);

        if(completeChains.isEmpty())
            return false;

        // take the shortest chain
        int shortestLength = -1;
        SvChain shortestChain = null;
        for(SvChain chain : completeChains)
        {
            if(chain.getLinkCount() < 2)
                continue; // only consider reciprocal chains if this short

            chain.recalcLength();

            if(shortestChain == null || chain.getLength() < shortestLength)
            {
                shortestChain = chain;
                shortestLength = chain.getLength();
            }
        }

        if(shortestChain != null)
        {
            shortestChain.setId(cluster.getChains().size());

            LOGGER.info("sample({}) cluster({}) adding complete chain({}) length({}) with {} linked pairs:",
                    sampleId, cluster.getId(), shortestChain.getId(), shortestChain.getLength(), shortestChain.getLinkCount());

            shortestChain.logLinks();

            cluster.addChain(shortestChain);
            return true;
        }

        return false;
    }

    private static boolean findCompletedLinks(
            final List<SvClusterData> varList, final List<SvLinkedPair> requiredLinks, final List<SvLinkedPair> existingLinks,
            final List<SvChain> partialChains, List<SvLinkedPair> testLinks, int currentIndex, List<SvChain> completeChains)
    {
        if(currentIndex >= testLinks.size())
            return false;

        int requiredSvCount = varList.size();

        LOGGER.debug("testLinkIndex({}) existingLinks({}) completedChains({}) partialChains({}) reqSVCount({})",
                currentIndex, existingLinks.size(), completeChains.size(), partialChains.size(), requiredSvCount);

        // try adding each of the remaining links recursively to get a set of completely-links breakends
        // and if so, add this collection of potential chains
        for(int i = currentIndex; i < testLinks.size(); ++i)
        {
            final SvLinkedPair testLink = testLinks.get(i);

            // start with existing partial chains and then build on them
            List<SvChain> workingChains = Lists.newArrayList();

            // must be a deep copy to avoid changing the set of working chains at this level of the search
            for(final SvChain chain : partialChains)
            {
                workingChains.add(new SvChain(chain));
            }

            List<SvLinkedPair> workingLinks = Lists.newArrayList();
            workingLinks.addAll(existingLinks);

            if(!canAddLinkedPair(workingLinks, workingChains, testLink, requiredSvCount))
                continue;

            // add this link to any chains if possible
            addLinkToChains(testLink, workingChains);

            // every new SV is now part of a partial chain (even if only a single link)
            workingLinks.add(testLink);

            // reconcile chains together if possible
            if(workingChains.size() > 1)
            {
                for (int firstIndex = 0; firstIndex < workingChains.size(); ++firstIndex)
                {
                    final SvChain firstChain = workingChains.get(firstIndex);

                    if(firstChain.getSvCount() == requiredSvCount || firstChain.isClosedLoop())
                        continue;

                    int nextIndex = firstIndex + 1;
                    while (nextIndex < workingChains.size())
                    {
                        final SvChain nextChain = workingChains.get(nextIndex);

                        if(nextChain.getSvCount() == requiredSvCount || nextChain.isClosedLoop())
                        {
                            ++nextIndex;
                            continue;
                        }

                        if (reconcileChains(firstChain, nextChain))
                        {
                            workingChains.remove(nextIndex);
                        }
                        else
                        {
                            ++nextIndex;
                        }
                    }
                }
            }

            // and check whether any chains are complete
            int chainIndex = 0;
            while(chainIndex < workingChains.size())
            {
                final SvChain chain = workingChains.get(chainIndex);

                if(chain.getSvCount() < requiredSvCount)
                {
                    ++chainIndex;
                    continue;
                }

                boolean isDuplicate = false;
                for(final SvChain completeChain : completeChains)
                {
                    if(completeChain.isIdentical(chain))
                    {
                        isDuplicate = true;

                        LOGGER.debug("skipping duplicate complete chain:");
                        chain.logLinks();
                        break;
                    }
                }

                boolean hasRequiredLinks = chain.hasLinks(requiredLinks);

                if(!isDuplicate && hasRequiredLinks)
                {
                    chain.setId(completeChains.size());
                    completeChains.add(chain);

                    LOGGER.debug("iters({} of {}) testLinkIndex({}) existingLinks({}) completedChains({}) workingChains({}) - adding complete chain:",
                            i, testLinks.size(), currentIndex, existingLinks.size(), completeChains.size(), workingChains.size());

                    // LOGGER.debug("added complete potential chain:");
                    chain.logLinks();
                }

                workingChains.remove(chainIndex);
            }

            if(workingChains.isEmpty())
            {
                // this search path is done so no need to keep searching for more links to add
                break;
            }

            LOGGER.debug("iters({} of {}) testLinkIndex({}) existingLinks({}) completedChains({}) workingChains({}) - continuing search",
                    i, testLinks.size(), currentIndex, existingLinks.size(), completeChains.size(), workingChains.size());

            // continue the search, moving on to try adding the next test link
            findCompletedLinks(varList, requiredLinks, workingLinks, workingChains, testLinks, i + 1, completeChains);
        }

        return true;
    }

    private static void addLinkToChains(SvLinkedPair linkedPair, List<SvChain> chains)
    {
        // add to an existing chain or create a new one
        for(final SvChain chain : chains)
        {
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

}
