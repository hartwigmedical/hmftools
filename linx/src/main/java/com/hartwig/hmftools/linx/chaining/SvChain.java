package com.hartwig.hmftools.linx.chaining;

import static java.lang.Math.abs;

import static com.hartwig.hmftools.common.utils.Strings.appendStr;
import static com.hartwig.hmftools.common.sv.StructuralVariantType.DUP;
import static com.hartwig.hmftools.linx.LinxConfig.LNX_LOGGER;
import static com.hartwig.hmftools.linx.analysis.SvUtilities.calcConsistency;
import static com.hartwig.hmftools.linx.chaining.ChainUtils.getSequenceStr;
import static com.hartwig.hmftools.common.utils.sv.StartEndIterator.SE_END;
import static com.hartwig.hmftools.common.utils.sv.StartEndIterator.SE_PAIR;
import static com.hartwig.hmftools.common.utils.sv.StartEndIterator.SE_START;
import static com.hartwig.hmftools.common.utils.sv.StartEndIterator.seIndex;
import static com.hartwig.hmftools.linx.chaining.LinkFinder.getMinTemplatedInsertionLength;
import static com.hartwig.hmftools.linx.chaining.LinkFinder.isPossibleLink;

import java.util.List;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.linx.types.SvBreakend;
import com.hartwig.hmftools.linx.types.LinkedPair;
import com.hartwig.hmftools.linx.types.SvVarData;

public class SvChain {

    private int mId;

    // links are added in such a way the the 'first' SV in the link is the link to the preceding
    // link in the chain, and the 'second' SV links to its other breakend in the next link in the chain
    private List<SvVarData> mSvList;
    private List<LinkedPair> mLinkedPairs;
    private SvBreakend[] mOpenBreakends; // cached for convenience
    private double mJcn;
    private double mJcnUncertainty;

    private boolean mIsClosedLoop;
    private boolean mDoubleMinute; // the forming one, not subclonal
    private int mLinkSum; // simple comparison check

    public SvChain(int chainId)
    {
        mId = chainId;
        mSvList = Lists.newArrayList();
        mLinkedPairs = Lists.newArrayList();
        mJcn = 0;
        mJcnUncertainty = 0;
        mIsClosedLoop = false;
        mDoubleMinute = false;
        mLinkSum = 0;
        mOpenBreakends = new SvBreakend[SE_PAIR];
    }

    public int id() { return mId; }
    public void setId(int id) { mId = id; }

    public int getSvCount() { return mSvList.size(); }
    public List<SvVarData> getSvList() { return mSvList; }

    public int getLinkCount() { return mLinkedPairs.size(); }
    public List<LinkedPair> getLinkedPairs() { return mLinkedPairs; }

    public void addLink(final LinkedPair pair, boolean addToStart)
    {
        if(!mSvList.contains(pair.first()))
            mSvList.add(pair.first());

        if(!mSvList.contains(pair.second()))
            mSvList.add(pair.second());

        LinkedPair newPair = LinkedPair.copy(pair);

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

    public boolean firstLinkOpenOnStart()
    {
        return !mLinkedPairs.isEmpty() ? mLinkedPairs.get(0).firstUnlinkedOnStart() : false;
    }
    public boolean lastLinkOpenOnStart() { return !mLinkedPairs.isEmpty() ? !mLinkedPairs.get(mLinkedPairs.size()-1).secondLinkOnStart() : false; }

    public String toString()
    {
        return String.format("%d: svs(%d) links(%d) jcn(%.1f) %s",
                mId, mSvList.size(), mLinkedPairs.size(), mJcn, mIsClosedLoop ? "closed" : "");
    }

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

    public boolean canAddLinkedPairToStart(final LinkedPair pair)
    {
        return canAddLinkedPair(pair, true, false);
    }

    public boolean canAddLinkedPairToEnd(final LinkedPair pair)
    {
        return canAddLinkedPair(pair, false, false);
    }

    public boolean canAddLinkedPair(final LinkedPair pair, boolean toStart, boolean allowReplicated)
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

    public boolean linkWouldCloseChain(final LinkedPair pair)
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

    public boolean couldCloseChain()
    {
        if(mSvList.size() == 1 && mSvList.get(0).type() == DUP)
            return true;

        final SvBreakend chainStart = getOpenBreakend(true);
        final SvBreakend chainEnd = getOpenBreakend(false);

        if(chainStart != null && !chainStart.getSV().isSglBreakend() && chainEnd != null && !chainEnd.getSV().isSglBreakend())
        {
            if(chainEnd.getSV() == chainStart.getSV() && chainStart != chainEnd)
            {
                return true;
            }
            else
            {
                int minTILength = getMinTemplatedInsertionLength(chainStart, chainEnd);

                return isPossibleLink(chainStart.chromosome(), chainStart.position(), chainStart.orientation(),
                        chainEnd.chromosome(), chainEnd.position(), chainEnd.orientation(), minTILength);
            }
        }

        return false;
    }

    public boolean isClosedLoop()
    {
        return mIsClosedLoop || (mSvList.size() == 1 && mSvList.get(0).type() == DUP);
    }

    public boolean isDoubleMinute() { return mDoubleMinute; }
    public void setDoubleMinute(boolean toggle) { mDoubleMinute = toggle; }

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

        LinkedPair link = LinkedPair.from(chainEnd, chainStart);
        link.setLinkReason(reason, linkIndex);

        mLinkedPairs.add(link);
        mIsClosedLoop = true;

        mLinkSum += chainEnd.getSV().id() + chainStart.getSV().id();
    }

    public void closeChain() { mIsClosedLoop = true; }

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

            for(final LinkedPair pair : mLinkedPairs)
            {
                if(pair.firstBreakend().arm() != pair.secondBreakend().arm())
                {
                    if(traversesCentromere)
                        return false;

                    traversesCentromere = true;
                }
            }

            if(traversesCentromere) // orientations don't matter if exactly one centromere is included
                return true;

            // must start and end on same arm with consistent breakends
            if(!mOpenBreakends[SE_START].getChrArm().equals(mOpenBreakends[SE_END].getChrArm()))
                return false;

            return mOpenBreakends[SE_START].orientation() != mOpenBreakends[SE_END].orientation();
        }
        else
        {
            return calcConsistency(mOpenBreakends[SE_START]) + calcConsistency(mOpenBreakends[SE_END]) == 0;
        }
    }


    public void copyFrom(final SvChain otherChain)
    {
        for(final LinkedPair pair : otherChain.getLinkedPairs())
        {
            addLink(pair, false);
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
            final LinkedPair pair = chain.getLinkedPairs().get(i);

            if(pair.first() == pair.second() && !pair.isDupLink())
                return false;

            if(!chain.getSvList().contains(pair.first()) || !chain.getSvList().contains(pair.second()))
                return false;

            final LinkedPair nextPair = chain.getLinkedPairs().get(i+1);

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
            final LinkedPair pair = mLinkedPairs.get(i);

            LNX_LOGGER.debug("chain({}) {}: pair({}) {} length({}) index({})",
                    mId, i, pair.toString(), pair.getLinkReason(), pair.baseLength(), pair.getLinkIndex());
        }
    }

    public long getLength(boolean closeEnds)
    {
        // defined as the sum of the TI lengths
        long length = mLinkedPairs.stream().mapToLong(x -> abs(x.baseLength())).sum();

        if(closeEnds)
        {
            final SvBreakend chainStart = getOpenBreakend(true);
            final SvBreakend chainEnd = getOpenBreakend(false);

            if(chainEnd != null && chainStart != null)
            {
                // skip the special single DUP case
                if(mLinkedPairs.size() == 1 && chainEnd.getSV() == chainStart.getSV())
                    return length;

                if(chainStart.chromosome().equals(chainEnd.chromosome()))
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
            final LinkedPair pair = mLinkedPairs.get(index);

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

    public int getAssemblyLinkCount() { return (int)mLinkedPairs.stream().filter(x -> x.isAssembled()).count(); }

    public boolean hasMatchingSVs(final SvChain other)
    {
        // return true if the 'other' chain's SVs are all in this chain (and it can have more)
        for(final SvVarData var : other.getSvList())
        {
            if(!mSvList.stream().anyMatch(x -> x == var))
                return false;
        }

        return true;
    }

    public boolean sameLinks(final SvChain other)
    {
        if(other.getLinkCount() != mLinkedPairs.size())
            return false;

        for(final LinkedPair pair : mLinkedPairs)
        {
            if(!other.getLinkedPairs().stream().anyMatch(x -> x.matches(pair)))
                return false;
        }

        return true;
    }

    public int linkSum() { return mLinkSum; }

}
