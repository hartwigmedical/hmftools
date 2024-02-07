package com.hartwig.hmftools.linx.chaining;

import static com.hartwig.hmftools.common.utils.sv.StartEndIterator.SE_END;
import static com.hartwig.hmftools.common.utils.sv.StartEndIterator.SE_START;
import static com.hartwig.hmftools.common.utils.sv.StartEndIterator.isStart;
import static com.hartwig.hmftools.linx.LinxConfig.LNX_LOGGER;
import static com.hartwig.hmftools.linx.analysis.SvUtilities.formatJcn;
import static com.hartwig.hmftools.linx.analysis.SvUtilities.makeChrArmStr;
import static com.hartwig.hmftools.linx.chaining.ChainJcnLimits.jcnMatch;

import java.util.List;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.linx.types.LinkedPair;
import com.hartwig.hmftools.linx.types.SvBreakend;
import com.hartwig.hmftools.linx.types.SvVarData;

public class ChainUtils
{
    public static final int CHAIN_LINK_COUNT = 0;
    public static final int CHAIN_ASSEMBLY_LINK_COUNT = 1;
    public static final int CHAIN_LENGTH = 2;

    public static int[] breakendsAreChained(final SvChain chain, final SvVarData var1, boolean v1Start, final SvVarData var2, boolean v2Start)
    {
        // check whether these breakends face towards each other in the chain
        // return info about the number of links including now many are assembled

        boolean be1FacesUp = false; // 'up' here means towards a higher index
        int be1Index = -1;
        boolean be2FacesUp = false;
        int be2Index = -1;

        int[] linkData = new int[CHAIN_LENGTH+1];

        // in every linked pair, the first element is lower in the chain (where 0 is considering the beginning) and the second is higher,
        // so for every SV it's 'first' breakend in a given pair faces up, the 'second' faces down
        for(int i = 0; i < chain.getLinkedPairs().size(); ++i)
        {
            final LinkedPair pair = chain.getLinkedPairs().get(i);

            if(pair.first() == var1 && pair.firstLinkOnStart() == v1Start)
            {
                be1Index = i;
                be1FacesUp = true;
            }
            else if(pair.second() == var1 && pair.secondLinkOnStart() == v1Start)
            {
                be1Index = i;
                be1FacesUp = false;
            }

            if(pair.first() == var2 && pair.firstLinkOnStart() == v2Start)
            {
                be2Index = i;
                be2FacesUp = true;
            }
            else if(pair.second() == var2 && pair.secondLinkOnStart() == v2Start)
            {
                be2Index = i;
                be2FacesUp = false;
            }

            if(be1Index >= 0 || be2Index >= 0)
            {
                ++linkData[CHAIN_LINK_COUNT];

                linkData[CHAIN_LENGTH] += pair.baseLength();

                if(pair.isAssembled())
                    ++linkData[CHAIN_ASSEMBLY_LINK_COUNT];
            }

            if(be1Index >=0 && be2Index >= 0)
            {
                if((be1Index <= be2Index && be1FacesUp && !be2FacesUp)
                        || (be1Index >= be2Index && !be1FacesUp && be2FacesUp))
                {
                    return linkData;
                }
            }
        }

        // reset since both links weren't found
        linkData[CHAIN_LINK_COUNT] = 0;
        linkData[CHAIN_ASSEMBLY_LINK_COUNT] = 0;

        return linkData;
    }

    public static void foldbackChainOnLink(final SvChain chain, final LinkedPair pair1, final LinkedPair pair2)
    {
        // validate the links being added
        final SvBreakend chainBreakend;
        if(pair1.first() == pair2.first() && pair1.second() == pair2.second())
        {
            chainBreakend = pair1.firstBreakend() == pair2.firstBreakend() ? pair1.firstBreakend() : pair1.secondBreakend();
        }
        else if(pair1.first() == pair2.second() && pair1.second() == pair2.first())
        {
            chainBreakend = pair1.firstBreakend() == pair2.secondBreakend() ? pair1.firstBreakend() : pair1.secondBreakend();
        }
        else
        {
            LNX_LOGGER.error("chain({}) failed to add foldback pairs: {} and {}", chain.id(), pair1, pair2);
            return;
        }

        boolean connectOnStart = chain.getOpenBreakend(true) == chainBreakend;
        int linkCount = chain.getLinkedPairs().size();

        if(connectOnStart)
        {
            List<LinkedPair> existingLinks = Lists.newArrayList(chain.getLinkedPairs());

            chain.addLink(pair1, true);
            chain.addLink(pair2, true);

            // the beginning of the chain will now form the middle
            for(LinkedPair pair : existingLinks)
            {
                chain.addLink(pair, true);
            }
        }
        else
        {
            chain.addLink(pair1, false);
            chain.addLink(pair2, false);

            // the beginning of the chain will now form the middle
            for(int index = linkCount - 1; index >= 0; --index)
            {
                final LinkedPair pair = chain.getLinkedPairs().get(index);
                chain.addLink(pair, false);
            }
        }
    }

    public static void foldbackChainOnChain(final SvChain chain, final SvChain foldbackChain, final LinkedPair pair1, final LinkedPair pair2)
    {
        // the provided chain splits and replicates this chain
        // keep existing links in place
        // add the first connecting link
        // then the foldback chain
        // then the next connecting link
        // then repeat this chain's links in reverse

        // establish what is connecting to what
        final SvBreakend chainStart = chain.getOpenBreakend(true);
        final SvBreakend chainEnd = chain.getOpenBreakend(false);

        final SvBreakend fbBreakendStart = foldbackChain.getOpenBreakend(true);
        final SvBreakend fbBreakendEnd = foldbackChain.getOpenBreakend(false);

        final SvBreakend otherBreakend = pair1.hasBreakend(fbBreakendStart) ?
                pair1.getOtherBreakend(fbBreakendStart) : pair1.getOtherBreakend(fbBreakendEnd);

        boolean connectOnStart = chainStart == otherBreakend;

        if((connectOnStart && !pair2.hasBreakend(chainStart)) || (!connectOnStart && !pair2.hasBreakend(chainEnd)))
        {
            LNX_LOGGER.error("chain({}) failed to add foldback pairs: {} and {}", chain.id(), pair1, pair2);
            return;
        }

        List<LinkedPair> existingLinks = Lists.newArrayList(chain.getLinkedPairs());

        final List<LinkedPair> fbLinks = foldbackChain.getLinkedPairs();

        // doesn't matter which pair is added first
        chain.addLink(pair1, connectOnStart);

        final SvBreakend connectingBreakend = connectOnStart ? chainStart : chainEnd;
        boolean connectingFoldbackChainStart = foldbackChain.getOpenBreakend(true) == pair1.getOtherBreakend(connectingBreakend);

        if(connectingFoldbackChainStart)
        {
            for(final LinkedPair fbPair : fbLinks)
            {
                chain.addLink(fbPair, connectOnStart);
            }
        }
        else
        {
            for(int index = fbLinks.size() - 1; index >= 0; --index)
            {
                final LinkedPair fbPair = fbLinks.get(index);
                chain.addLink(fbPair, connectOnStart);
            }
        }

        // now add the second pair
        chain.addLink(pair2, connectOnStart);

        // and then repeat this chain's links in reverse from before
        // the beginning of the chain will now form the middle
        if(!connectOnStart)
        {
            for(int index = existingLinks.size() - 1; index >= 0; --index)
            {
                final LinkedPair pair = existingLinks.get(index);
                chain.addLink(pair, connectOnStart);
            }
        }
        else
        {
            for(final LinkedPair pair : existingLinks)
            {
                chain.addLink(pair, connectOnStart);
            }
        }
    }

    public static void duplicateChainOnLink(final SvChain chain, final LinkedPair pair1, final LinkedPair pair2)
    {
        // validate the links being added
        final LinkedPair endPair;
        final LinkedPair startPair;
        if(chain.canAddLinkedPairToEnd(pair1) && chain.canAddLinkedPairToStart(pair2))
        {
            endPair = pair1;
            startPair = pair2;
        }
        else if(chain.canAddLinkedPairToEnd(pair2) && chain.canAddLinkedPairToStart(pair1))
        {
            endPair = pair2;
            startPair = pair1;
        }
        else
        {
            return;
        }

        int linkCount = chain.getLinkedPairs().size();

        chain.addLink(endPair, false);
        chain.addLink(startPair, false);

        // the beginning of the chain will now form the middle
        for(int index = 0; index < linkCount; ++index)
        {
            final LinkedPair pair = chain.getLinkedPairs().get(index);
            chain.addLink(pair, false);
        }
    }

    public static boolean reverseSectionOnBreakend(final SvChain chain, final SvBreakend breakend)
    {
        // reverses a section of a chain
        // eg sAe - B - eCs - D if reversed at C's start breakend will become sCe - B - eAs - D

        for(int i = 0; i < chain.getLinkedPairs().size(); ++i)
        {
            final LinkedPair pair = chain.getLinkedPairs().get(i);

            // find the location in the chain where the breakend is joined and join its closest chain end here instead
            if(pair.firstBreakend() == breakend)
            {
                SvBreakend chainStart = chain.getOpenBreakend(true);
                LinkedPair newLink = LinkedPair.from(chainStart, pair.secondBreakend());

                if(newLink.positionDistance() <= 0)
                    return false;

                newLink.setLinkReason("RECIP_INV_RECONFIG", 0);

                List<LinkedPair> linksToSwitch = Lists.newArrayList();
                for(int j = 0; j < i; ++j)
                {
                    linksToSwitch.add(chain.getLinkedPairs().get(j));
                }

                for(int j = 0; j < linksToSwitch.size() + 1; ++j)
                {
                    chain.getLinkedPairs().remove(0);
                }

                chain.addLink(newLink, true);

                for(int j = 0; j < linksToSwitch.size(); ++j)
                {
                    chain.addLink(linksToSwitch.get(j), true);
                }

                break;
            }
            else if(pair.secondBreakend() == breakend)
            {
                SvBreakend chainEnd = chain.getOpenBreakend(false);
                LinkedPair newLink = LinkedPair.from(pair.firstBreakend(), chainEnd);

                if(newLink.positionDistance() <= 0)
                    return false;

                newLink.setLinkReason("RECIP_INV_RECONFIG", 0);

                List<LinkedPair> linksToSwitch = Lists.newArrayList();
                for(int j = i + 1; j < chain.getLinkedPairs().size(); ++j)
                {
                    linksToSwitch.add(chain.getLinkedPairs().get(j));
                }

                for(int j = 0; j < linksToSwitch.size() + 1; ++j)
                {
                    chain.getLinkedPairs().remove(i);
                }

                chain.addLink(newLink, false);

                for(int j = linksToSwitch.size() - 1; j >= 0; --j)
                {
                    chain.addLink(linksToSwitch.get(j), false);
                }

                break;
            }
        }

        return true;
    }

    public static void reconcileChains(final List<SvChain> chains)
    {
        reconcileChains(chains, false, 0, false);
    }

    public static void reconcileChains(final List<SvChain> chains, boolean checkChainSplits, int nextChainId, boolean useChainEndJCNs)
    {
        // join 2 chains into a single chain if they have the same SV's opposing breakends on their ends
        // if the chains differ in JCN then either skip merging them or split off (copy) the larger JCN chain before joining
        if(chains.size() <= 1)
            return;

        // look for chains with opposite breakends of the same SV
        int index1 = 0;
        boolean skippedChainClosing = false;

        while(index1 < chains.size() - 1)
        {
            SvChain chain1 = chains.get(index1);

            boolean chainsMerged = false;

            for(int index2 = index1 + 1; index2 < chains.size(); ++index2)
            {
                SvChain chain2 = chains.get(index2);

                if(chain1 == chain2)
                    continue;

                if(checkChainSplits && identicalChain(chain1, chain2, false, false))
                {
                    chain1.setJcnData(chain1.jcn() + chain2.jcn(), chain1.jcnUncertainty());
                    chains.remove(index2);
                    chainsMerged = true;
                    break;
                }

                // avoid merging chains which close a loop
                /*
                if(chain2.getOpenBreakend(true) != null && chain2.getOpenBreakend(false) != null)
                {
                    final SvLinkedPair chain2Pair = SvLinkedPair.from(
                            chain2.getOpenBreakend(true).getOtherBreakend(), chain2.getOpenBreakend(false).getOtherBreakend());

                    if(chain1.linkWouldCloseChain(chain2Pair))
                    {
                        skippedChainClosing = true;
                        continue;
                    }
                }
                */

                boolean jcnMatched = jcnMatch(chain1.jcn(), chain1.jcnUncertainty(), chain2.jcn(), chain2.jcnUncertainty());

                if(!jcnMatched && !checkChainSplits && !useChainEndJCNs)
                    continue;

                for(int be1 = SE_START; be1 <= SE_END; ++be1)
                {
                    boolean c1Start = isStart(be1);

                    final SvBreakend breakend1 = chain1.getOpenBreakend(c1Start);

                    if(breakend1 == null)
                        continue;

                    for(int be2 = SE_START; be2 <= SE_END; ++be2)
                    {
                        boolean c2Start = isStart(be2);

                        final SvBreakend breakend2 = chain2.getOpenBreakend(c2Start);

                        if(breakend2 == null)
                            continue;

                        boolean couldJoinChains = breakend1 != breakend2 && breakend1.getSV() == breakend2.getSV();

                        if(!couldJoinChains)
                            continue;

                        if(useChainEndJCNs && !jcnMatched)
                        {
                            // even if the chains don't have a match, check the breakends themselves vs the chains
                            if(jcnMatch(chain1.jcn(), chain1.jcnUncertainty(), breakend1.jcn(), breakend1.jcnUncertainty())
                                    && jcnMatch(chain2.jcn(), chain2.jcnUncertainty(), breakend1.jcn(), breakend1.jcnUncertainty()))
                            {
                                jcnMatched = true;
                            }
                            else if(!checkChainSplits)
                            {
                                continue;
                            }
                        }

                        if(!jcnMatched)
                        {
                            // lower the JCN of the higher chain and add its copy to the end of the chains list
                            SvChain newChain = new SvChain(nextChainId++);

                            SvChain higherJcnChain = chain1.jcn() > chain2.jcn() ? chain1 : chain2;
                            SvChain lowerJcnChain = higherJcnChain == chain1 ? chain2 : chain1;

                            newChain.copyFrom(higherJcnChain);
                            newChain.setJcnData(higherJcnChain.jcn() - lowerJcnChain.jcn(), higherJcnChain.jcnUncertainty());

                            LNX_LOGGER.debug("splitting chain({}) jcn({}) vs chain({}) jcn({}) into new chain({}) JCN({})",
                                    higherJcnChain.id(), formatJcn(higherJcnChain.jcn()),
                                    lowerJcnChain.id(), formatJcn(lowerJcnChain.jcn()),
                                    newChain.id(), formatJcn(newChain.jcn()));

                            higherJcnChain.setJcnData(lowerJcnChain.jcn(), higherJcnChain.jcnUncertainty());

                            chains.add(newChain);
                        }

                        LNX_LOGGER.debug("merging chain({} links={}) {} to chain({} links={}) {}",
                                chain1.id(), chain1.getLinkCount(), c1Start ? "start" : "end",
                                chain2.id(), chain2.getLinkCount(), c2Start ? "start" : "end");

                        if(c2Start)
                        {
                            // merge chains and remove the latter
                            for(LinkedPair linkedPair : chain2.getLinkedPairs())
                            {
                                chain1.addLink(linkedPair, c1Start);
                            }
                        }
                        else
                        {
                            // add in reverse
                            for(int index = chain2.getLinkedPairs().size() - 1; index >= 0; --index)
                            {
                                LinkedPair linkedPair = chain2.getLinkedPairs().get(index);
                                chain1.addLink(linkedPair, c1Start);
                            }
                        }

                        chains.remove(index2);

                        chainsMerged = true;
                        break;
                    }

                    if(chainsMerged)
                        break;
                }

                if(chainsMerged)
                    break;
            }

            if(!chainsMerged)
            {
                ++index1;
            }
            else if(skippedChainClosing)
            {
                // need to retry from the beginning
                skippedChainClosing = false;
                index1 = 0;

            }
        }
    }

    public static boolean identicalChain(final SvChain chain, final SvChain other, boolean allowSubsets)
    {
        return identicalChain(chain, other, allowSubsets, false);
    }

    public static boolean identicalChain(final SvChain chain, final SvChain otherChain, boolean allowSubsets, boolean allowSameLinks)
    {
        if(chain.linkSum() != otherChain.linkSum())
            return false;

        /*
        // form a chain in reverse to compare

        final SvChain otherChain = other;
        if(!reversed)
        {
            otherChain = other;
        }
        else
        {
            otherChain = new SvChain(other.id());
            for(int i = mLinkedPairs.size() - 1; i >= 0; --i)
            {
                final SvLinkedPair pair = other.getLinkedPairs().get(i);
                SvLinkedPair reversePair = SvLinkedPair.from(pair.secondBreakend(), pair.firstBreakend());
                otherChain.addLink(reversePair, false);
            }
        }
        */

        // same SVs forming same links in the same order
        final List<LinkedPair> otherLinks = otherChain.getLinkedPairs();

        if(chain.getLinkedPairs().size() == otherLinks.size())
        {
            boolean exactMatch = true;

            for(int i = 0; i < chain.getLinkedPairs().size(); ++i)
            {
                final LinkedPair pair = chain.getLinkedPairs().get(i);
                final LinkedPair otherPair = otherLinks.get(i);

                if(pair.first() != otherPair.first() || pair.second() != otherPair.second())
                {
                    exactMatch = false;
                    break;
                }
            }

            if(exactMatch)
                return true;

            return allowSameLinks && chain.sameLinks(otherChain);
        }
        else
        {
            if(!allowSubsets)
                return false;

            final List<LinkedPair> longerLinks = chain.getLinkedPairs().size() > otherLinks.size() ? chain.getLinkedPairs() : otherLinks;
            final List<LinkedPair> shorterLinks = chain.getLinkedPairs().size() < otherLinks.size() ? chain.getLinkedPairs() : otherLinks;

            int j = 0;
            boolean matchingLinkFound = false;
            for(int i = 0; i < longerLinks.size(); ++i)
            {
                final LinkedPair pair = longerLinks.get(i);

                if(j >= shorterLinks.size())
                    break;

                final LinkedPair otherPair = shorterLinks.get(j);

                boolean linksMatch = pair.first() == otherPair.first() && pair.second() == otherPair.second();

                if(linksMatch)
                {
                    matchingLinkFound = true;
                    ++j;
                }
                else if(matchingLinkFound)
                {
                    // links matched up to this point but not differ
                    return false;
                }
            }

            return matchingLinkFound;
        }
    }


    private static final String CHAIN_SEQ_DELIM = " - ";

    private static String breakendSeqStr(final SvBreakend breakend)
    {
        return String.format("%s_%s_%s",
                breakend.usesStart() ? "s" : "e", breakend.getSV().id(), breakend.usesStart() ? "e" : "s");
    }

    public static String getSequenceStr(final SvChain chain)
    {
        String sequenceStr = "";

        for(int i = 0; i < chain.getLinkedPairs().size(); ++i)
        {
            LinkedPair pair = chain.getLinkedPairs().get(i);

            if(i == 0)
            {
                SvBreakend startBreakend = chain.getOpenBreakend(true);

                if(startBreakend != null)
                {
                    sequenceStr += makeChrArmStr(startBreakend.chromosome(), startBreakend.arm()) + "_" + startBreakend.direction();
                    sequenceStr += CHAIN_SEQ_DELIM;

                    sequenceStr += breakendSeqStr(startBreakend);
                    sequenceStr += CHAIN_SEQ_DELIM;
                }
                else
                {
                    sequenceStr += "sgl_unclear";

                    sequenceStr += CHAIN_SEQ_DELIM;

                    sequenceStr += breakendSeqStr(pair.first().getBreakend(true));
                    sequenceStr += CHAIN_SEQ_DELIM;
                }
            }

            sequenceStr += breakendSeqStr(pair.secondBreakend());
            sequenceStr += CHAIN_SEQ_DELIM;

            if(i == chain.getLinkedPairs().size() - 1)
            {
                SvBreakend endBreakend = chain.getOpenBreakend(false);
                if(endBreakend != null)
                {
                    sequenceStr += makeChrArmStr(endBreakend.chromosome(), endBreakend.arm()) + "_" + endBreakend.direction();
                }
                else
                {
                    sequenceStr += "sgl_unclear";
                }
            }
        }

        return sequenceStr;
    }

}
