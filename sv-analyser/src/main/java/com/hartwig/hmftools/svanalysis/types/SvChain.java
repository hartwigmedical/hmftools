package com.hartwig.hmftools.svanalysis.types;

import static java.lang.Math.abs;

import static com.hartwig.hmftools.svanalysis.analysis.LinkFinder.areLinkedSection;
import static com.hartwig.hmftools.svanalysis.analysis.SvUtilities.appendStr;
import static com.hartwig.hmftools.svanalysis.analysis.SvUtilities.calcConsistency;
import static com.hartwig.hmftools.svanalysis.analysis.SvUtilities.makeChrArmStr;

import java.util.List;

import com.google.common.collect.Lists;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;

public class SvChain {

    private int mId;

    // links are added in such a way the the 'first' SV in the link is the link to the preceding
    // link in the chain, and the 'second' SV links to its other breakend in the next link in the chain
    private List<SvVarData> mSvList;
    private List<SvLinkedPair> mLinkedPairs;

    private boolean mIsClosedLoop;
    private boolean mIsValid;
    private int mReplicationCount; // if a chain is identically repeated, report the number of times this occurs

    private String mDetails;

    private static final Logger LOGGER = LogManager.getLogger(SvChain.class);

    public SvChain(int chainId)
    {
        mId = chainId;
        mSvList = Lists.newArrayList();
        mLinkedPairs = Lists.newArrayList();
        mIsClosedLoop = false;
        mIsValid = true;
        mDetails = "";
        mReplicationCount = 1;
    }

    public int id() { return mId; }
    public void setId(int id) { mId = id; }

    public int getSvCount() { return mSvList.size(); }
    public List<SvVarData> getSvList() { return mSvList; }

    public int getLinkCount() { return mLinkedPairs.size(); }
    public List<SvLinkedPair> getLinkedPairs() { return mLinkedPairs; }

    public void addLink(final SvLinkedPair pair, boolean addToStart)
    {
        if(mLinkedPairs.isEmpty())
        {
            mLinkedPairs.add(pair);

            mSvList.add(pair.first());
            mSvList.add(pair.second());
            return;
        }

        if(mLinkedPairs.contains(pair))
            return;

        // check ordering and switch if required so that the 'first' SV always links to the preceding link and vice versa
        if((addToStart && pair.second() != mLinkedPairs.get(0).first())
        || (!addToStart && pair.first() != mLinkedPairs.get(mLinkedPairs.size()-1).second()))
        {
            pair.switchSVs();
        }

        if(addToStart)
            mLinkedPairs.add(0, pair); // insert at front
        else
            mLinkedPairs.add(pair);

        final SvVarData first = pair.first();
        final SvVarData second = pair.second();

        boolean containsFirst = mSvList.contains(first);
        boolean containsSecond = mSvList.contains(second);

        int lastIndex = mSvList.size() - 1;

        if(containsFirst && containsSecond)
        {
            // check that these SVs are at the start and end, otherwise the new link is invalid
            if((mSvList.get(0) == first && mSvList.get(lastIndex) != second) || (mSvList.get(0) == second && mSvList.get(lastIndex) != first))
            {
                mIsValid = false;
                return;
            }

            // no need to add an SV twice (ie to both start and end)
            mIsClosedLoop = true;
        }
        else
        {
            if (addToStart)
            {
                if (mSvList.get(0) == first || mSvList.get(1) == first)
                {
                    // second SV is the new one here
                    mSvList.add(0, second);
                }
                else
                {
                    mSvList.add(0, first);
               }
            }
            else
            {
                if (mSvList.get(lastIndex-1) == first || mSvList.get(lastIndex) == first)
                {
                    mSvList.add(second);
                }
                else
                {
                    mSvList.add(first);
                }
            }
        }

        setIsValid();
    }

    public SvVarData getFirstSV() { return mSvList.isEmpty() ? null : mSvList.get(0); }
    public SvVarData getLastSV() { return mSvList.isEmpty() ? null : mSvList.get(mSvList.size()-1); }
    public SvVarData getChainEndSV(boolean isFirst) { return isFirst ? getFirstSV() : getLastSV(); }

    public SvLinkedPair getLinkedPair(boolean isStart) { return isStart ? getFirstLinkedPair() : getLastLinkedPair(); }
    public SvLinkedPair getFirstLinkedPair() { return mLinkedPairs.isEmpty() ? null : mLinkedPairs.get(0); }
    public SvLinkedPair getLastLinkedPair() { return mLinkedPairs.isEmpty() ? null : mLinkedPairs.get(mLinkedPairs.size()-1); }

    public boolean firstLinkOpenOnStart()
    {
        return !mLinkedPairs.isEmpty() ? mLinkedPairs.get(0).firstUnlinkedOnStart() : false;
    }
    public boolean lastLinkOpenOnStart() { return !mLinkedPairs.isEmpty() ? !mLinkedPairs.get(mLinkedPairs.size()-1).secondLinkOnStart() : false; }
    public boolean chainEndOpenOnStart(boolean isFirst) { return isFirst ? firstLinkOpenOnStart() : lastLinkOpenOnStart(); }

    public final SvBreakend getOpenBreakend(boolean isStart)
    {
        return isStart ? getFirstSV().getBreakend(firstLinkOpenOnStart()) : getLastSV().getBreakend(lastLinkOpenOnStart());
    }

    public boolean isClosedLoop() { return mIsClosedLoop; }

    public boolean hasReplicatedSVs()
    {
        for(int i = 0; i < mSvList.size(); ++i)
        {
            final SvVarData var1 = mSvList.get(i);

            for(int j = i+1; j < mSvList.size(); ++j)
            {
                if (mSvList.get(j).equals(var1, true))
                    return true;
            }
        }

        return false;
    }

    public boolean isConsistent()
    {
        if(getFirstSV().isNullBreakend() || getLastSV().isNullBreakend())
            return false;

        // treat the start and end breakends like those of a single SV
        int consistency = calcConsistency(getFirstSV(), firstLinkOpenOnStart());
        consistency += calcConsistency(getLastSV(), lastLinkOpenOnStart());
        return consistency == 0;
    }

    public String getDetails() { return mDetails; }
    public void setDetails(final String details) { mDetails = details; }

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

        if(pair.first().equals(chainSv, allowReplicated) && pair.firstUnlinkedOnStart() != chainOpenSide)
            return true;
        else if(pair.second().equals(chainSv, allowReplicated) && pair.secondUnlinkedOnStart() != chainOpenSide)
            return true;
        else
            return false;
    }

    public boolean linkWouldCloseChain(final SvLinkedPair pair)
    {
        return canAddLinkedPairToStart(pair) && canAddLinkedPairToEnd(pair);
    }

    public void setIsValid()
    {
        if(mSvList.isEmpty() || mLinkedPairs.isEmpty())
        {
            mIsValid = false;
            return;
        }

        // check for duplicate breakends
        for(int i = 0; i < mLinkedPairs.size(); ++i)
        {
            final SvLinkedPair lp1 = mLinkedPairs.get(i);

            for(int j = i + 1; j< mLinkedPairs.size(); ++j)
            {
                final SvLinkedPair lp2 = mLinkedPairs.get(j);

                if(lp1.hasLinkClash(lp2))
                {
                    mIsValid = false;
                    return;
                }
            }
        }

        mIsValid = true;
    }

    public int getReplicationCount() { return mReplicationCount; }
    public void addToReplicationCount() { ++mReplicationCount; }

    public void logLinks()
    {
        if(mLinkedPairs.size() < 50)
        {
            LOGGER.debug("chain({}): {}", mId, getSequenceStr(this));
        }

        for(int i = 0; i < mLinkedPairs.size(); ++i)
        {
            final SvLinkedPair pair = mLinkedPairs.get(i);

            LOGGER.debug("chain({}) {}: pair({}) {} {} length({})",
                    mId, i, pair.toString(), pair.assemblyInferredStr(), pair.getLinkReason(), pair.length());
        }
    }

    public int getLength(boolean closeEnds)
    {
        // defined as the sum of the TI lengths
        int length = 0;

        for(final SvLinkedPair pair : mLinkedPairs)
        {
            length += abs(pair.length());
        }

        if(closeEnds)
        {
            final SvBreakend chainStart = getOpenBreakend(true);
            final SvBreakend chainEnd = getOpenBreakend(false);

            if(mLinkedPairs.size() == 1 && chainEnd.getSV() == chainStart.getSV()) // skip the special single DUP case
                return length;

            if (chainStart.chromosome().equals(chainEnd.chromosome()))
            {
                length += abs(chainStart.position() - chainEnd.position());
            }
        }

        return length;
    }

    public int getSvCount(boolean includeReplicated)
    {
        if(includeReplicated)
            return mSvList.size();
        else
            return (int)mSvList.stream().filter(x -> !x.isReplicatedSv()).count();
    }

    public boolean hasSV(final SvVarData var, boolean allowReplicated)
    {
        for(final SvVarData chainSv : mSvList)
        {
            if (chainSv.equals(var, allowReplicated))
                return true;
        }

        return false;
    }

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

            if(pair.first().equals(var,true))
            {
                linkInfo = String.format("%d%s", index, pair.firstLinkOnStart() ? "s" : "e");
            }
            else if(pair.second().equals(var,true))
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

            if(pair.first().equals(var1, true) && pair.firstLinkOnStart() == v1Start)
            {
                be1Index = i;
                be1FacesUp = true;
            }
            else if(pair.second().equals(var1, true) && pair.secondLinkOnStart() == v1Start)
            {
                be1Index = i;
                be1FacesUp = false;
            }

            if(pair.first().equals(var2, true) && pair.firstLinkOnStart() == v2Start)
            {
                be2Index = i;
                be2FacesUp = true;
            }
            else if(pair.second().equals(var2, true) && pair.secondLinkOnStart() == v2Start)
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

    public static void checkChainReplication(final SvCluster cluster)
    {
        if(!cluster.hasReplicatedSVs())
            return;

        // check whether any chains are identical to others using replicated SVs
        // in which case remove the replicated SVs and chains
        List<SvChain> chains = cluster.getChains();

        int index1 = 0;
        while(index1 < chains.size())
        {
            final SvChain chain1 = chains.get(index1);

            for(int index2 = index1+1; index2 < chains.size(); ++index2)
            {
                final SvChain chain2 = chains.get(index2);

                if(chain1.identicalChain(chain2))
                {
                    boolean allReplicatedSVs = chain2.getSvCount() == 0;

                    LOGGER.debug("cluster({}) removing duplicate chain({}) vs origChain({}) all replicated({})",
                            cluster.id(), chain2.id(), chain1.id(), allReplicatedSVs);

                    // remove these replicated SVs as well as the replicated chain
                    if(allReplicatedSVs)
                    {
                        for (final SvVarData var : chain2.getSvList())
                        {
                            cluster.removeReplicatedSv(var);
                        }
                    }

                    chains.remove(index2);
                    continue;

                }

                ++index2;
            }

            ++index1;
        }
    }

    public boolean identicalChain(final SvChain other)
    {
        // same SVs forming same links in the same order
        if(mLinkedPairs.size() != other.getLinkCount())
            return false;

        for(int i = 0; i < mLinkedPairs.size(); ++i)
        {
            final SvLinkedPair pair = mLinkedPairs.get(i);
            final SvLinkedPair otherPair = other.getLinkedPairs().get(i);

            if(!pair.first().equals(otherPair.first(), true))
                return false;

            if(!pair.second().equals(otherPair.second(), true))
                return false;
        }

        return true;
    }

    private static String CHAIN_SEQ_DELIM = " - ";

    private static String breakendSeqStr(final SvBreakend breakend)
    {
        return String.format("%s_%s_%s",
                breakend.usesStart() ? "s" : "e", breakend.getOrigSV().id(), breakend.usesStart() ? "e" : "s");
    }

    public static String getSequenceStr(final SvChain chain)
    {
        String sequenceStr = "";

        for(int i = 0; i < chain.getLinkedPairs().size(); ++i)
        {
            SvLinkedPair pair = chain.getLinkedPairs().get(i);

            if(i == 0)
            {
                if(!pair.first().isNullBreakend())
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

            sequenceStr += breakendSeqStr(pair.getSecondBreakend());
            sequenceStr += CHAIN_SEQ_DELIM;

            if(i == chain.getLinkedPairs().size() - 1)
            {
                if(!pair.second().isNullBreakend())
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

    public static List<SvVarData> getRepeatedSvSequence(final List<SvVarData> svList, int firstIndex, int secondIndex, boolean walkForwards)
    {
        // walk forward from these 2 start points comparing SVs
        List<SvVarData> sequence = Lists.newArrayList();

        int i = firstIndex;
        int j = secondIndex + 1;

        if(walkForwards)
            ++i;
        else
            --i;

        while(i < secondIndex && i >= 0 && j < svList.size())
        {
            final SvVarData var1 = svList.get(i);
            final SvVarData var2 = svList.get(j);

            if(!var1.equals(var2, true))
                break;

            sequence.add(var1);

            ++j;

            if(walkForwards)
                ++i;
            else
                --i;
        }

        return sequence;
    }
}
