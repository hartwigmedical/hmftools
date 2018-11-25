package com.hartwig.hmftools.svanalysis.types;

import static java.lang.Math.abs;

import static com.hartwig.hmftools.svanalysis.analysis.SvUtilities.calcConsistency;

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

    // has an entry for each SV, indicating which end of the SV links to the preceding SV
    private List<Boolean> mSvLinkToPrecedingOnStart;

    private int mLength;
    private boolean mIsClosedLoop;
    private boolean mIsValid;

    private int mConsistencyCount;

    private static final Logger LOGGER = LogManager.getLogger(SvChain.class);

    public SvChain(int chainId)
    {
        mId = chainId;
        mSvList = Lists.newArrayList();
        mLinkedPairs = Lists.newArrayList();
        mSvLinkToPrecedingOnStart = Lists.newArrayList();
        mLength = 0;
        mIsClosedLoop = false;
        mIsValid = true;
        mConsistencyCount = 0;
    }

    public SvChain(final SvChain other)
    {
        mId = other.getId();
        mSvList = Lists.newArrayList();
        mLinkedPairs = Lists.newArrayList();
        mSvLinkToPrecedingOnStart = Lists.newArrayList();
        mLength = 0;
        mIsClosedLoop = false;
        mIsValid = true;

        for(final SvLinkedPair pair : other.getLinkedPairs())
        {
            addLink(pair, canAddLinkedPairToStart(pair));
        }
    }

    public int getId() { return mId; }
    public void setId(int id) { mId = id; }

    public boolean isValid() { return mIsValid; }

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
            mSvLinkToPrecedingOnStart.add(!pair.firstLinkOnStart());
            mSvList.add(pair.second());

            // if this second variant in the linked pair is connected on its Start breakend, then
            // the precending SV link is 'Start'
            // if the other end of this variant is added as a new link, it will be on the other breakend,
            // so for now the end of the chain is open on the opposite of this last boolean value
            mSvLinkToPrecedingOnStart.add(pair.secondLinkOnStart());
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

                    // whichever end is not linked in this pair is the one exposed at the start of the chain
                    mSvLinkToPrecedingOnStart.add(0, !pair.secondLinkOnStart());
                }
                else
                {
                    mSvList.add(0, first);
                    mSvLinkToPrecedingOnStart.add(0, !pair.firstLinkOnStart());
                }
            }
            else
            {
                if (mSvList.get(lastIndex-1) == first || mSvList.get(lastIndex) == first)
                {
                    mSvList.add(second);
                    mSvLinkToPrecedingOnStart.add(pair.secondLinkOnStart());
                }
                else
                {
                    mSvList.add(first);
                    mSvLinkToPrecedingOnStart.add(pair.firstLinkOnStart());
                }
            }
        }

        setIsValid();

        mConsistencyCount = calcConsistency(mSvList);
    }

    public SvVarData getFirstSV() { return mSvList.isEmpty() ? null : mSvList.get(0); }
    public SvVarData getLastSV() { return mSvList.isEmpty() ? null : mSvList.get(mSvList.size()-1); }

    public SvLinkedPair getFirstLinkedPair() { return mLinkedPairs.isEmpty() ? null : mLinkedPairs.get(0); }
    public SvLinkedPair getLastLinkedPair() { return mLinkedPairs.isEmpty() ? null : mLinkedPairs.get(mLinkedPairs.size()-1); }

    public boolean firstLinkOpenOnStart()
    {
        return !mLinkedPairs.isEmpty() ? mLinkedPairs.get(0).firstUnlinkedOnStart() : false;
    }

    public boolean lastLinkOpenOnStart()
    {
        return !mLinkedPairs.isEmpty() ? !mLinkedPairs.get(mLinkedPairs.size()-1).secondLinkOnStart() : false;
    }

    /*
    public boolean firstLinkOpenOnStart() { return mSvLinkToPrecedingOnStart.isEmpty() ? false : mSvLinkToPrecedingOnStart.get(0); }

    public boolean lastLinkOpenOnStart()
    {
        if(mSvLinkToPrecedingOnStart.isEmpty())
            return false;

        // the array of booleans refers to the the link to the preceding SV,
        // so if it's linked on the start backwards, the end must be open and vice versa
        return !mSvLinkToPrecedingOnStart.get(mSvList.size()-1);
    }
    */

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

    public boolean isConsistent() { return mConsistencyCount == 0; }

    public boolean canAddLinkedPairToStart(final SvLinkedPair pair)
    {
        if(mLinkedPairs.isEmpty())
            return true;

        if(mLinkedPairs.contains(pair))
            return false;

        // check for the same SV exposed in opposite ends in the 2 linked pairs
        if(pair.first() == getFirstSV() && pair.firstUnlinkedOnStart() != firstLinkOpenOnStart())
            return true;
        else if(pair.second() == getFirstSV() && pair.secondUnlinkedOnStart() != firstLinkOpenOnStart())
            return true;
        else
            return false;
    }

    public boolean canAddLinkedPairToEnd(final SvLinkedPair pair)
    {
        if(mLinkedPairs.isEmpty())
            return false;

        if(mLinkedPairs.contains(pair))
            return false;

        if(pair.first() == getLastSV() && pair.firstUnlinkedOnStart() != lastLinkOpenOnStart())
            return true;
        else if(pair.second() == getLastSV() && pair.secondUnlinkedOnStart() != lastLinkOpenOnStart())
            return true;
        else
            return false;
    }

    public boolean linkWouldCloseChain(final SvLinkedPair pair)
    {
        return canAddLinkedPairToStart(pair) && canAddLinkedPairToEnd(pair);
    }

    public boolean isIdentical(final SvChain other)
    {
        // true if both closed loops and containing the same linked pairs
        // or if not closed, having the same end points
        if(mLinkedPairs.size() != other.getLinkedPairs().size())
            return false;

        List<SvLinkedPair> otherLinkedPairs = Lists.newArrayList();
        otherLinkedPairs.addAll(other.getLinkedPairs());

        for(final SvLinkedPair pair : mLinkedPairs)
        {
            boolean matched = false;
            for(int i = 0; i < otherLinkedPairs.size(); ++i)
            {
                final SvLinkedPair otherPair = otherLinkedPairs.get(i);
                if(pair.matches(otherPair, true))
                {
                    otherLinkedPairs.remove(i); // reduce each time to make subsequent matches faster
                    matched = true;
                    break;
                }
            }

            if(!matched)
                return false;
        }

        if(isClosedLoop() && other.isClosedLoop())
            return true;

        if(!getFirstSV().equals(other.getFirstSV(), true) || firstLinkOpenOnStart() != other.firstLinkOpenOnStart())
            return false;

        if(!getLastSV().equals(other.getLastSV(), true) || lastLinkOpenOnStart() != other.lastLinkOpenOnStart())
            return false;

        return true;
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

    public void logLinks()
    {
        for(int i = 0; i < mLinkedPairs.size(); ++i)
        {
            final SvLinkedPair pair = mLinkedPairs.get(i);

            LOGGER.debug("chain({}) {}: pair({}) {} {} len={}",
                    mId, i, pair.toString(), pair.linkType(), pair.isInferred() ? "inferred" : "assembly", pair.length());
        }
    }

    public int getLength() { return mLength; }

    public void recalcLength()
    {
        // defined as the TI and DB lengths plus the variant lengths for those linked at both end
        mLength = 0;

        for(int i = 0; i < mLinkedPairs.size(); ++i)
        {
            final SvLinkedPair pair = mLinkedPairs.get(i);

            mLength += abs(pair.length());

            if(i > 0 && i < mLinkedPairs.size())
            {
                mLength += pair.first().length();
            }
        }
    }

    public int getUniqueSvCount()
    {
        int count = 0;

        for(final SvVarData var : mSvList)
        {
            if(!var.isReplicatedSv())
                ++count;
        }

        return count;
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
        for(int index = 0; index < mSvList.size(); ++index)
        {
            if(!mSvList.get(index).equals(var, true))
                continue;

            if(!varIndices.isEmpty())
                varIndices += ";";

            varIndices += index;
        }

        return varIndices;
    }

    public int getSvIndex(final SvVarData var, boolean matchStart)
    {
        for(int index = 0; index < mLinkedPairs.size(); ++index)
        {
            final SvLinkedPair pair = mLinkedPairs.get(index);

            if((pair.first().equals(var) && pair.firstLinkOnStart() != matchStart)
            || (pair.second().equals(var) && pair.secondLinkOnStart() != matchStart))
            {
                return index;
            }
        }

        return -1;
    }

    public boolean breakendsAreChained(final SvVarData var1, boolean v1Start, final SvVarData var2, boolean v2Start)
    {
        // the 2 boolean passed in are the ends which are linked through the chain

        // check whether these breakends face towards each other in the chain
        boolean be1FacesUp = false; // 'up' here means towards a higher index
        int be1Index = -1;
        boolean be2FacesUp = false;
        int be2Index = -1;

        for(int i = 0; i < mSvList.size(); ++i)
        {
            final SvVarData var = mSvList.get(i);

            if(var.equals(var1, true))
            {
                be1Index = i;

                // if start links to preceding then end is facing up, so true if is checking the end
                be1FacesUp = mSvLinkToPrecedingOnStart.get(i) != v1Start;
            }
            else if(var.equals(var2, true))
            {
                be2Index = i;
                be2FacesUp = mSvLinkToPrecedingOnStart.get(i) != v2Start;
            }

            if(be1Index >=0 && be2Index >= 0)
            {
                if(be1Index < be2Index && be1FacesUp && !be2FacesUp)
                    return true;
                else if(be1Index > be2Index && !be1FacesUp && be2FacesUp)
                    return true;
            }
        }

        return false;
    }

    public boolean hasLinkedPair(final SvLinkedPair pair)
    {
        return mLinkedPairs.contains(pair);
    }

}
