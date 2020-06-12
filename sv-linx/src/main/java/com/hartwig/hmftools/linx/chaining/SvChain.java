package com.hartwig.hmftools.linx.chaining;

import static java.lang.Math.abs;

import static com.hartwig.hmftools.common.utils.Strings.appendStr;
import static com.hartwig.hmftools.common.variant.structural.StructuralVariantType.DUP;
import static com.hartwig.hmftools.common.variant.structural.StructuralVariantType.SGL;
import static com.hartwig.hmftools.linx.LinxConfig.LNX_LOGGER;
import static com.hartwig.hmftools.linx.analysis.SvUtilities.calcConsistency;
import static com.hartwig.hmftools.linx.analysis.SvUtilities.formatJcn;
import static com.hartwig.hmftools.linx.analysis.SvUtilities.makeChrArmStr;
import static com.hartwig.hmftools.linx.chaining.ChainJcnLimits.jcnMatch;
import static com.hartwig.hmftools.linx.types.SvLinkedPair.LOCATION_TYPE_EXTERNAL;
import static com.hartwig.hmftools.linx.types.SvLinkedPair.LOCATION_TYPE_INTERNAL;
import static com.hartwig.hmftools.linx.types.SvLinkedPair.LOCATION_TYPE_REMOTE;
import static com.hartwig.hmftools.common.utils.sv.StartEndIterator.SE_END;
import static com.hartwig.hmftools.common.utils.sv.StartEndIterator.SE_PAIR;
import static com.hartwig.hmftools.common.utils.sv.StartEndIterator.SE_START;
import static com.hartwig.hmftools.common.utils.sv.StartEndIterator.isStart;
import static com.hartwig.hmftools.common.utils.sv.StartEndIterator.seIndex;
import static com.hartwig.hmftools.linx.types.SvConstants.SHORT_TI_LENGTH;

import java.util.List;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.linx.types.SvBreakend;
import com.hartwig.hmftools.linx.types.SvLinkedPair;
import com.hartwig.hmftools.linx.types.SvVarData;

public class SvChain {

    private int mId;

    // links are added in such a way the the 'first' SV in the link is the link to the preceding
    // link in the chain, and the 'second' SV links to its other breakend in the next link in the chain
    private List<SvVarData> mSvList;
    private List<SvLinkedPair> mLinkedPairs;
    private SvBreakend[] mOpenBreakends; // cached for convenience
    private double mJcn;
    private double mJcnUncertainty;

    private boolean mIsClosedLoop;
    private int mLinkSum; // simple comparison check

    public SvChain(int chainId)
    {
        mId = chainId;
        mSvList = Lists.newArrayList();
        mLinkedPairs = Lists.newArrayList();
        mJcn = 0;
        mJcnUncertainty = 0;
        mIsClosedLoop = false;
        mLinkSum = 0;
        mOpenBreakends = new SvBreakend[SE_PAIR];
    }

    public int id() { return mId; }
    public void setId(int id) { mId = id; }

    public int getSvCount() { return mSvList.size(); }
    public List<SvVarData> getSvList() { return mSvList; }

    public int getLinkCount() { return mLinkedPairs.size(); }
    public List<SvLinkedPair> getLinkedPairs() { return mLinkedPairs; }

    public void addLink(final SvLinkedPair pair, boolean addToStart)
    {
        if(!mSvList.contains(pair.first()))
            mSvList.add(pair.first());

        if(!mSvList.contains(pair.second()))
            mSvList.add(pair.second());

        SvLinkedPair newPair = SvLinkedPair.copy(pair);

        mLinkSum += pair.first().id() + pair.second().id();

        if(mLinkedPairs.isEmpty())
        {
            mLinkedPairs.add(newPair);
            updateBreakends();
            return;
        }

        // validate can be added
        if((addToStart && !canAddLinkedPairToStart(newPair)) || (!addToStart && !canAddLinkedPairToEnd(newPair)))
        {
            LNX_LOGGER.error("chain({}) trying to add invalid link({}) to {}", mId, newPair.toString(), addToStart ? "start" : "end");
            return;
        }

        // check ordering and switch if required so that the 'first' SV always links to the preceding link and vice versa
        if((addToStart && newPair.second() != mLinkedPairs.get(0).first())
        || (!addToStart && newPair.first() != mLinkedPairs.get(mLinkedPairs.size() - 1).second()))
        {
            newPair.switchSVs();
        }

        if(addToStart)
            mLinkedPairs.add(0, newPair); // insert at front
        else
            mLinkedPairs.add(newPair);

        updateBreakends();
    }

    public void setJcnData(double jcn, double uncertainty)
    {
        mJcn = jcn;
        mJcnUncertainty = uncertainty;
    }

    public double jcn() { return mJcn; }
    public double jcnUncertainty() { return mJcnUncertainty; }

    public SvVarData getChainEndSV(boolean isStart)
    {
        if(mLinkedPairs.isEmpty())
            return null;

        if(isStart)
            return mLinkedPairs.get(0).first();
        else
            return mLinkedPairs.get(mLinkedPairs.size() - 1).second();
    }

    public SvVarData getFirstSV() { return getChainEndSV(true); }
    public SvVarData getLastSV() { return getChainEndSV(false); }

    public SvLinkedPair getLinkedPair(boolean isStart)
    {
        if(mLinkedPairs.isEmpty())
            return null;

        if(isStart)
            return mLinkedPairs.get(0);
        else
            return mLinkedPairs.get(mLinkedPairs.size()-1);
    }

    public boolean firstLinkOpenOnStart()
    {
        return !mLinkedPairs.isEmpty() ? mLinkedPairs.get(0).firstUnlinkedOnStart() : false;
    }
    public boolean lastLinkOpenOnStart() { return !mLinkedPairs.isEmpty() ? !mLinkedPairs.get(mLinkedPairs.size()-1).secondLinkOnStart() : false; }

    private void updateBreakends()
    {
        mOpenBreakends[SE_START] = getFirstSV().getBreakend(firstLinkOpenOnStart());
        mOpenBreakends[SE_END] = getLastSV().getBreakend(lastLinkOpenOnStart());
    }

    public final SvBreakend getOpenBreakend(int seIndex) { return mOpenBreakends[seIndex]; }
    public final SvBreakend getOpenBreakend(boolean isStart)
    {
        return getOpenBreakend(seIndex(isStart));
    }

    public boolean isConsistent() { return isConsistent(false); }

    public boolean isConsistent(boolean checkCentromere)
    {
        if(getFirstSV().isSglBreakend() || getLastSV().isSglBreakend())
            return false;

        if(checkCentromere)
        {
            // criteria for consistency - either starts and ends on same arm with consistent breakends OR
            // if ends on a different arm it needs to go through a centromere once and have 2 telomeres
            boolean traversesCentromere = false;

            for (final SvLinkedPair pair : mLinkedPairs)
            {
                if (pair.firstBreakend().arm() != pair.secondBreakend().arm())
                {
                    if (traversesCentromere)
                        return false;

                    traversesCentromere = true;
                }
            }

            if (traversesCentromere) // orientations don't matter if exactly one centromere is included
                return true;

            // must start and end on same arm with consistent breakends
            if (!mOpenBreakends[SE_START].getChrArm().equals(mOpenBreakends[SE_END].getChrArm()))
                return false;

            return mOpenBreakends[SE_START].orientation() != mOpenBreakends[SE_END].orientation();
        }
        else
        {
            return calcConsistency(mOpenBreakends[SE_START]) + calcConsistency(mOpenBreakends[SE_END]) == 0;
        }


    }

    public boolean canAddLinkedPairToStart(final SvLinkedPair pair)
    {
        return canAddLinkedPair(pair, true, false);
    }

    public boolean canAddLinkedPairToEnd(final SvLinkedPair pair)
    {
        return canAddLinkedPair(pair, false, false);
    }

    public boolean canAddLinkedPair(final SvLinkedPair pair, boolean toStart, boolean allowReplicated)
    {
        if(mLinkedPairs.isEmpty())
            return true;

        if(mLinkedPairs.contains(pair))
            return false;

        // check for the same SV exposed in opposite ends in the 2 linked pairs
        final SvVarData chainSv = toStart ? getFirstSV() : getLastSV();
        boolean chainOpenSide = toStart ? firstLinkOpenOnStart() : lastLinkOpenOnStart();

        if(pair.first() == chainSv && pair.firstUnlinkedOnStart() != chainOpenSide)
            return true;
        else if(pair.second() == chainSv && pair.secondUnlinkedOnStart() != chainOpenSide)
            return true;
        else
            return false;
    }

    public boolean linkWouldCloseChain(final SvLinkedPair pair)
    {
        final SvBreakend chainStart = getOpenBreakend(true);
        final SvBreakend chainEnd = getOpenBreakend(false);

        if(pair.firstBreakend() == chainStart && pair.secondBreakend() == chainEnd)
            return true;
        else if(pair.firstBreakend() == chainEnd && pair.secondBreakend() == chainStart)
            return true;
        else
            return false;
    }

    public boolean isClosedLoop()
    {
        return mIsClosedLoop || (mSvList.size() == 1 && mSvList.get(0).type() == DUP);
    }

    public void closeChain(final String reason, int linkIndex)
    {
        final SvBreakend chainStart = getOpenBreakend(true);
        final SvBreakend chainEnd = getOpenBreakend(false);

        if(!chainStart.chromosome().equals(chainEnd.chromosome()))
        {
            LNX_LOGGER.error("chain({}) cannot close chain with ends: {} & {}",
                    mId, chainStart.toString(), chainEnd.toString());
            return;
        }

        SvLinkedPair link = SvLinkedPair.from(chainEnd, chainStart);
        link.setLinkReason(reason, linkIndex);

        mLinkedPairs.add(link);
        mIsClosedLoop = true;

        mLinkSum += chainEnd.getSV().id() + chainStart.getSV().id();
    }

    public void closeChain() { mIsClosedLoop = true; }

    public void foldbackChainOnLink(final SvLinkedPair pair1, final SvLinkedPair pair2)
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
            LNX_LOGGER.error("chain({}) failed to add foldback pairs: {} and {}", mId, pair1, pair2);
            return;
        }

        boolean connectOnStart = getOpenBreakend(true) == chainBreakend;
        int linkCount = mLinkedPairs.size();

        if(connectOnStart)
        {
            List<SvLinkedPair> existingLinks = Lists.newArrayList(mLinkedPairs);

            addLink(pair1, true);
            addLink(pair2, true);

            // the beginning of the chain will now form the middle
            for(SvLinkedPair pair : existingLinks)
            {
                addLink(pair, true);
            }
        }
        else
        {
            addLink(pair1, false);
            addLink(pair2, false);

            // the beginning of the chain will now form the middle
            for(int index = linkCount - 1; index >= 0; --index)
            {
                final SvLinkedPair pair = mLinkedPairs.get(index);
                addLink(pair, false);
            }
        }
    }

    public void foldbackChainOnChain(final SvChain foldbackChain, final SvLinkedPair pair1, final SvLinkedPair pair2)
    {
        // the provided chain splits and replicates this chain
        // keep existing links in place
        // add the first connecting link
        // then the foldback chain
        // then the next connecting link
        // then repeat this chain's links in reverse

        // establish what is connecting to what
        final SvBreakend chainStart = getOpenBreakend(true);
        final SvBreakend chainEnd = getOpenBreakend(false);

        final SvBreakend fbBreakendStart = foldbackChain.getOpenBreakend(true);
        final SvBreakend fbBreakendEnd = foldbackChain.getOpenBreakend(false);

        final SvBreakend otherBreakend = pair1.hasBreakend(fbBreakendStart) ?
                pair1.getOtherBreakend(fbBreakendStart) : pair1.getOtherBreakend(fbBreakendEnd);

        boolean connectOnStart = chainStart == otherBreakend;

        if((connectOnStart && !pair2.hasBreakend(chainStart)) || (!connectOnStart && !pair2.hasBreakend(chainEnd)))
        {
            LNX_LOGGER.error("chain({}) failed to add foldback pairs: {} and {}", mId, pair1, pair2);
            return;
        }

        List<SvLinkedPair> existingLinks = Lists.newArrayList(mLinkedPairs);

        final List<SvLinkedPair> fbLinks = foldbackChain.getLinkedPairs();

        // doesn't matter which pair is added first
        addLink(pair1, connectOnStart);

        final SvBreakend connectingBreakend = connectOnStart ? chainStart : chainEnd;
        boolean connectingFoldbackChainStart = foldbackChain.getOpenBreakend(true) == pair1.getOtherBreakend(connectingBreakend);

        if(connectingFoldbackChainStart)
        {
            for (final SvLinkedPair fbPair : fbLinks)
            {
                addLink(fbPair, connectOnStart);
            }
        }
        else
        {
            for(int index = fbLinks.size() - 1; index >= 0; --index)
            {
                final SvLinkedPair fbPair = fbLinks.get(index);
                addLink(fbPair, connectOnStart);
            }
        }

        // now add the second pair
        addLink(pair2, connectOnStart);

        // and then repeat this chain's links in reverse from before
        // the beginning of the chain will now form the middle
        if(!connectOnStart)
        {
            for (int index = existingLinks.size() - 1; index >= 0; --index)
            {
                final SvLinkedPair pair = existingLinks.get(index);
                addLink(pair, connectOnStart);
            }
        }
        else
        {
            for(final SvLinkedPair pair : existingLinks)
            {
                addLink(pair, connectOnStart);
            }
        }
    }

    public void copyFrom(final SvChain otherChain)
    {
        for(final SvLinkedPair pair : otherChain.getLinkedPairs())
        {
            addLink(pair, false);
        }
    }

    public void duplicateChainOnLink(final SvLinkedPair pair1, final SvLinkedPair pair2)
    {
        // validate the links being added
        final SvLinkedPair endPair;
        final SvLinkedPair startPair;
        if(canAddLinkedPairToEnd(pair1) && canAddLinkedPairToStart(pair2))
        {
            endPair = pair1;
            startPair = pair2;
        }
        else if(canAddLinkedPairToEnd(pair2) && canAddLinkedPairToStart(pair1))
        {
            endPair = pair2;
            startPair = pair1;
        }
        else
        {
            return;
        }

        int linkCount = mLinkedPairs.size();

        addLink(endPair, false);
        addLink(startPair, false);

        // the beginning of the chain will now form the middle
        for(int index = 0; index < linkCount; ++index)
        {
            final SvLinkedPair pair = mLinkedPairs.get(index);
            addLink(pair, false);
        }
    }

    public boolean reverseSectionOnBreakend(final SvBreakend breakend)
    {
        // reverses a section of a chain
        // eg sAe - B - eCs - D if reversed at C's start breaknd will become sCe - B - eAs - D

        for(int i = 0; i < mLinkedPairs.size(); ++i)
        {
            final SvLinkedPair pair = mLinkedPairs.get(i);

            // find the location in the chain where the breakend is joined and join its closest chain end here instead
            if(pair.firstBreakend() == breakend)
            {
                SvBreakend chainStart = getOpenBreakend(true);
                SvLinkedPair newLink = SvLinkedPair.from(chainStart, pair.secondBreakend());

                if(newLink.length() <= 0)
                    return false;

                newLink.setLinkReason("RECIP_INV_RECONFIG", 0);

                List<SvLinkedPair> linksToSwitch = Lists.newArrayList();
                for(int j = 0; j < i; ++j)
                {
                    linksToSwitch.add(mLinkedPairs.get(j));
                }

                for(int j = 0; j < linksToSwitch.size() + 1; ++j)
                {
                    mLinkedPairs.remove(0);
                }

                addLink(newLink, true);

                for(int j = 0; j < linksToSwitch.size(); ++j)
                {
                    addLink(linksToSwitch.get(j), true);
                }

                break;
            }
            else if(pair.secondBreakend() == breakend)
            {
                SvBreakend chainEnd = getOpenBreakend(false);
                SvLinkedPair newLink = SvLinkedPair.from(pair.firstBreakend(), chainEnd);

                if(newLink.length() <= 0)
                    return false;

                newLink.setLinkReason("RECIP_INV_RECONFIG", 0);

                List<SvLinkedPair> linksToSwitch = Lists.newArrayList();
                for(int j = i + 1; j < mLinkedPairs.size(); ++j)
                {
                    linksToSwitch.add(mLinkedPairs.get(j));
                }

                for(int j = 0; j < linksToSwitch.size() + 1; ++j)
                {
                    mLinkedPairs.remove(i);
                }

                addLink(newLink, false);

                for(int j = linksToSwitch.size() - 1; j >= 0; --j)
                {
                    addLink(linksToSwitch.get(j), false);
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

    public static void reconcileChains(final List<SvChain> chains, boolean checkChainSplits, int nextChainId, boolean useChainEndPloidies)
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

            for (int index2 = index1 + 1; index2 < chains.size(); ++index2)
            {
                SvChain chain2 = chains.get(index2);

                if(chain1 == chain2)
                    continue;

                if(checkChainSplits && chain1.identicalChain(chain2, false, false))
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

                    if (chain1.linkWouldCloseChain(chain2Pair))
                    {
                        skippedChainClosing = true;
                        continue;
                    }
                }
                */

                boolean jcnMatched = jcnMatch(chain1.jcn(), chain1.jcnUncertainty(), chain2.jcn(), chain2.jcnUncertainty());

                if(!jcnMatched && !checkChainSplits && !useChainEndPloidies)
                    continue;

                for (int be1 = SE_START; be1 <= SE_END; ++be1)
                {
                    boolean c1Start = isStart(be1);

                    final SvBreakend breakend1 = chain1.getOpenBreakend(c1Start);

                    if(breakend1 == null)
                        continue;

                    for (int be2 = SE_START; be2 <= SE_END; ++be2)
                    {
                        boolean c2Start = isStart(be2);

                        final SvBreakend breakend2 = chain2.getOpenBreakend(c2Start);

                        if(breakend2 == null)
                            continue;

                        boolean couldJoinChains = breakend1 != breakend2 && breakend1.getSV() == breakend2.getSV();

                        if(!couldJoinChains)
                            continue;

                        if(useChainEndPloidies && !jcnMatched)
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

                        chains.remove(index2);

                        chainsMerged = true;
                        break;
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
            else if(skippedChainClosing)
            {
                // need to retry from the beginning
                skippedChainClosing = false;
                index1 = 0;

            }
        }
    }

    public static boolean checkIsValid(final SvChain chain)
    {
        if(chain.getSvCount() == 0 || chain.getLinkCount() == 0)
        {
            return false;
        }

        // check for duplicate breakends
        for(int i = 0; i < chain.getLinkCount() - 1; ++i)
        {
            final SvLinkedPair pair = chain.getLinkedPairs().get(i);

            if(pair.first() == pair.second() && !pair.isDupLink())
                return false;

            if(!chain.getSvList().contains(pair.first()) || !chain.getSvList().contains(pair.second()))
                return false;

            final SvLinkedPair nextPair = chain.getLinkedPairs().get(i+1);

            if(pair.second() != nextPair.first())
                return false;

            if(pair.secondLinkOnStart() == nextPair.firstLinkOnStart())
                return false;
        }

        return true;
    }

    public void logLinks()
    {
        if(mLinkedPairs.size() < 50)
        {
            LNX_LOGGER.debug("chain({}): {}", mId, getSequenceStr(this));
        }

        for(int i = 0; i < mLinkedPairs.size(); ++i)
        {
            final SvLinkedPair pair = mLinkedPairs.get(i);

            LNX_LOGGER.debug("chain({}) {}: pair({}) {} length({}) index({})",
                    mId, i, pair.toString(), pair.getLinkReason(), pair.length(), pair.getLinkIndex());
        }
    }

    public long getLength(boolean closeEnds)
    {
        // defined as the sum of the TI lengths
        long length = mLinkedPairs.stream().mapToLong(x -> abs(x.length())).sum();

        if(closeEnds)
        {
            final SvBreakend chainStart = getOpenBreakend(true);
            final SvBreakend chainEnd = getOpenBreakend(false);

            if(chainEnd != null && chainStart != null)
            {
                // skip the special single DUP case
                if (mLinkedPairs.size() == 1 && chainEnd.getSV() == chainStart.getSV())
                    return length;

                if (chainStart.chromosome().equals(chainEnd.chromosome()))
                {
                    length += abs(chainStart.position() - chainEnd.position());
                }
            }
        }

        return length;
    }

    public boolean hasSV(final SvVarData var) { return mSvList.contains(var); }

    public int getSvIndex(final SvVarData var)
    {
        for(int index = 0; index < mSvList.size(); ++index)
        {
            if(mSvList.get(index) == var)
                return index;
        }

        return -1;
    }

    public final String getSvIndices(final SvVarData var)
    {
        String varIndices = "";

        for(int index = 0; index < mLinkedPairs.size(); ++index)
        {
            final SvLinkedPair pair = mLinkedPairs.get(index);

            String linkInfo = "";

            if(pair.first() == var)
            {
                linkInfo = String.format("%d%s", index, pair.firstLinkOnStart() ? "s" : "e");
            }
            else if(pair.second() == var)
            {
                linkInfo = String.format("%d%s", index, pair.secondLinkOnStart() ? "s" : "e");
            }
            else
            {
                continue;
            }

            varIndices = appendStr(varIndices, linkInfo, ';');
        }

        return varIndices;
    }

    public boolean hasRepeatedSV()
    {
        if(mIsClosedLoop)
            return mLinkedPairs.size() > mSvList.size();
        else
            return mLinkedPairs.size() >= mSvList.size();
    }

    public int getAssemblyLinkCount()
    {
        int count = 0;
        for(final SvLinkedPair pair : mLinkedPairs)
        {
            if(pair.isAssembled())
                ++count;
        }

        return count;
    }

    public static int CHAIN_LINK_COUNT = 0;
    public static int CHAIN_ASSEMBLY_LINK_COUNT = 1;
    public static int CHAIN_LENGTH = 2;

    public int[] breakendsAreChained(final SvVarData var1, boolean v1Start, final SvVarData var2, boolean v2Start)
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
        for(int i = 0; i < mLinkedPairs.size(); ++i)
        {
            final SvLinkedPair pair = mLinkedPairs.get(i);

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

                linkData[CHAIN_LENGTH] += pair.length();

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

    public boolean identicalChain(final SvChain other, boolean allowSubsets)
    {
        return identicalChain(other, allowSubsets, false);
    }

    public boolean sameLinks(final SvChain other)
    {
        if(other.getLinkCount() != mLinkedPairs.size())
            return false;

        for(final SvLinkedPair pair : mLinkedPairs)
        {
            if(!other.getLinkedPairs().stream().anyMatch(x -> x.matches(pair)))
                return false;
        }

        return true;
    }

    public int linkSum() { return mLinkSum; }

    public boolean identicalChain(final SvChain otherChain, boolean allowSubsets, boolean allowSameLinks)
    {
        if(mLinkSum != otherChain.linkSum())
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
        final List<SvLinkedPair> otherLinks = otherChain.getLinkedPairs();

        if(mLinkedPairs.size() == otherLinks.size())
        {
            boolean exactMatch = true;

            for(int i = 0; i < mLinkedPairs.size(); ++i)
            {
                final SvLinkedPair pair = mLinkedPairs.get(i);
                final SvLinkedPair otherPair = otherLinks.get(i);

                if(pair.first() != otherPair.first() || pair.second() != otherPair.second())
                {
                    exactMatch = false;
                    break;
                }
            }

            if(exactMatch)
                return true;

            return allowSameLinks && sameLinks(otherChain);
        }
        else
        {
            if(!allowSubsets)
                return false;

            final List<SvLinkedPair> longerLinks = mLinkedPairs.size() > otherLinks.size() ? mLinkedPairs : otherLinks;
            final List<SvLinkedPair> shorterLinks = mLinkedPairs.size() < otherLinks.size() ? mLinkedPairs : otherLinks;

            int j = 0;
            boolean matchingLinkFound = false;
            for(int i = 0; i < longerLinks.size(); ++i)
            {
                final SvLinkedPair pair = longerLinks.get(i);

                if(j >= shorterLinks.size())
                    break;

                final SvLinkedPair otherPair = shorterLinks.get(j);

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

    private static String CHAIN_SEQ_DELIM = " - ";

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
            SvLinkedPair pair = chain.getLinkedPairs().get(i);

            if(i == 0)
            {
                if(!pair.first().isSglBreakend())
                {
                    SvBreakend startBreakend = chain.getOpenBreakend(true);

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
                if(!pair.second().isSglBreakend())
                {
                    SvBreakend endBreakend = chain.getOpenBreakend(false);
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

    public ChainMetrics extractChainMetrics()
    {
        final SvBreakend chainStart = getOpenBreakend(true);
        final SvBreakend chainEnd = getOpenBreakend(false);

        ChainMetrics metrics = new ChainMetrics();

        if(chainStart == null || chainEnd == null)
            return metrics;

        final SvBreakend lowerBreakend = chainStart.position() < chainEnd.position() ? chainStart : chainEnd;
        final SvBreakend upperBreakend = chainStart == lowerBreakend ? chainEnd : chainStart;

        boolean startEndSameArm = chainEnd.getChrArm().equals(chainStart.getChrArm());

        if(startEndSameArm)
        {
            if (lowerBreakend.orientation() == 1 && upperBreakend.orientation() == -1)
            {
                ++metrics.ChainEndsAway;
            }
            else if (lowerBreakend.orientation() == -1 && upperBreakend.orientation() == 1)
            {
                ++metrics.ChainEndsFace;
            }
        }

        for(final SvLinkedPair pair : mLinkedPairs)
        {
            if(pair.first().type() == SGL || pair.second().type() == SGL)
                continue;

            if(pair.locationType() == LOCATION_TYPE_INTERNAL)
            {
                ++metrics.InternalTIs;

                if(pair.length() <= SHORT_TI_LENGTH)
                    ++metrics.InternalShortTIs;

                if(pair.hasCopyNumberGain())
                    ++metrics.InternalTICnGain;
            }
            else if(pair.locationType() == LOCATION_TYPE_REMOTE || pair.locationType() == LOCATION_TYPE_EXTERNAL)
            {
                ++metrics.ExternalTIs;

                if(pair.length() <= SHORT_TI_LENGTH)
                    ++metrics.ExternalShortTIs;

                if(pair.hasCopyNumberGain())
                    ++metrics.ExternalTICnGain;
            }

            if(pair.overlapCount() > 0)
                ++metrics.OverlappingTIs;
        }

        return metrics;
    }

}
