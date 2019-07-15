package com.hartwig.hmftools.linx.chaining;

import static java.lang.Math.abs;
import static java.lang.Math.max;

import static com.hartwig.hmftools.common.variant.structural.StructuralVariantType.DUP;
import static com.hartwig.hmftools.common.variant.structural.StructuralVariantType.SGL;
import static com.hartwig.hmftools.linx.analysis.SvUtilities.appendStr;
import static com.hartwig.hmftools.linx.analysis.SvUtilities.calcConsistency;
import static com.hartwig.hmftools.linx.analysis.SvUtilities.copyNumbersEqual;
import static com.hartwig.hmftools.linx.analysis.SvUtilities.makeChrArmStr;
import static com.hartwig.hmftools.linx.types.SvLinkedPair.LOCATION_TYPE_INTERNAL;
import static com.hartwig.hmftools.linx.types.SvLinkedPair.LOCATION_TYPE_REMOTE;
import static com.hartwig.hmftools.linx.types.SvaConstants.SHORT_TI_LENGTH;

import java.util.List;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.linx.types.SvBreakend;
import com.hartwig.hmftools.linx.types.SvLinkedPair;
import com.hartwig.hmftools.linx.types.SvVarData;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;

public class SvChain {

    private int mId;

    // links are added in such a way the the 'first' SV in the link is the link to the preceding
    // link in the chain, and the 'second' SV links to its other breakend in the next link in the chain
    private List<SvVarData> mSvList;
    private List<SvLinkedPair> mLinkedPairs;
    private double mPloidy;
    private double mPloidyUncertainty;

    private boolean mIsClosedLoop;

    private static final Logger LOGGER = LogManager.getLogger(SvChain.class);

    public SvChain(int chainId)
    {
        mId = chainId;
        mSvList = Lists.newArrayList();
        mLinkedPairs = Lists.newArrayList();
        mPloidy = 0;
        mPloidyUncertainty = 0;
        mIsClosedLoop = false;
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

            if(!mSvList.contains(pair.second()))
                mSvList.add(pair.second());

            return;
        }

        if(mLinkedPairs.contains(pair))
            return;

        // check ordering and switch if required so that the 'first' SV always links to the preceding link and vice versa
        if((addToStart && pair.second() != mLinkedPairs.get(0).first())
        || (!addToStart && pair.first() != mLinkedPairs.get(mLinkedPairs.size() - 1).second()))
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
        boolean sameSV = (first == second);

        if(mSvList.isEmpty())
        {
            LOGGER.error("invalid empty list");
            return;
        }

        // it's possible for the list to only contain a single SV (if it's a DM DUP)
        int lastIndex = mSvList.size() - 1;
        int secondLastIndex = mSvList.size() >= 2 ? mSvList.size() - 2 : 0;
        int secondIndex = mSvList.size() >= 2 ? 1 : 0;

        if(containsFirst && containsSecond && !sameSV)
        {
            // check that these SVs are at the start and end, otherwise the new link is invalid
            if((mSvList.get(0) == first && mSvList.get(lastIndex) != second) || (mSvList.get(0) == second && mSvList.get(lastIndex) != first))
            {
                LOGGER.error("cannot add new pair: {}", pair.toString());
                return;
            }

            // no need to add an SV twice (ie to both start and end)
            mIsClosedLoop = true;
        }
        else
        {
            if (addToStart)
            {
                if (mSvList.get(0) == first || mSvList.get(secondIndex) == first)
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
                if (mSvList.get(secondLastIndex) == first || mSvList.get(lastIndex) == first)
                {
                    mSvList.add(second);
                }
                else
                {
                    mSvList.add(first);
                }
            }
        }
    }

    public void addLink(final SvLinkedPair pair, int index)
    {
        if(index >= mLinkedPairs.size())
            return;

        mLinkedPairs.add(index, pair);
    }

    public void setPloidyData(double ploidy, double uncertainty)
    {
        mPloidy = ploidy;
        mPloidyUncertainty = uncertainty;
    }

    public double ploidy() { return mPloidy; }
    public double ploidyUncertainty() { return mPloidyUncertainty; }

    public SvVarData getChainEndSV(boolean isFirst)
    {
        if(mSvList.isEmpty())
            return null;

        if (isFirst)
            return mSvList.get(0);
        else
            return mSvList.get(mSvList.size() - 1);
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

    public final SvBreakend getOpenBreakend(boolean isStart)
    {
        return isStart ? getFirstSV().getBreakend(firstLinkOpenOnStart()) : getLastSV().getBreakend(lastLinkOpenOnStart());
    }

    public boolean isClosedLoop()
    {
        return mIsClosedLoop || (mSvList.size() == 1 && mSvList.get(0).type() == DUP);
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

    public void foldbackChainOnLink(final SvLinkedPair pair1, final SvLinkedPair pair2)
    {
        // validate the links being added
        final SvVarData chainVar;
        if(pair1.first() == pair2.first() && pair1.second() == pair2.second())
        {
            chainVar = pair1.getFirstBreakend() == pair2.getFirstBreakend() ? pair1.first(): pair1.second();
        }
        else if(pair1.first() == pair2.second() && pair1.second() == pair2.first())
        {
            chainVar = pair1.getFirstBreakend() == pair2.getSecondBreakend() ? pair1.first(): pair1.second();
        }
        else
        {
            return;
        }

        boolean connectOnStart = chainVar == getFirstSV();
        int linkCount = mLinkedPairs.size();

        if(connectOnStart)
        {
            addLink(pair1, true);
            addLink(pair2, true);

            // the beginning of the chain will now form the middle
            for(int index = 0; index < linkCount; ++index)
            {
                final SvLinkedPair pair = mLinkedPairs.get(2 + index * 2);
                addLink(SvLinkedPair.from(pair.getFirstBreakend(), pair.getSecondBreakend()), true);
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
                addLink(SvLinkedPair.from(pair.getFirstBreakend(), pair.getSecondBreakend()), false);
            }
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
            addLink(SvLinkedPair.from(pair.getFirstBreakend(), pair.getSecondBreakend()), false);
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
            LOGGER.debug("chain({}): {}", mId, getSequenceStr(this));
        }

        for(int i = 0; i < mLinkedPairs.size(); ++i)
        {
            final SvLinkedPair pair = mLinkedPairs.get(i);

            LOGGER.debug("chain({}) {}: pair({}) {} {} length({}) index({})",
                    mId, i, pair.toString(), pair.assemblyInferredStr(), pair.getLinkReason(), pair.length(), pair.getLinkIndex());
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

    public boolean identicalChain(final SvChain other, boolean allowSubsets)
    {
        // same SVs forming same links in the same order
        final List<SvLinkedPair> otherLinks = other.getLinkedPairs();

        if(mLinkedPairs.size() == otherLinks.size())
        {
            for(int i = 0; i < mLinkedPairs.size(); ++i)
            {
                final SvLinkedPair pair = mLinkedPairs.get(i);
                final SvLinkedPair otherPair = otherLinks.get(i);

                if(!pair.first().equals(otherPair.first(), true))
                    return false;

                if(!pair.second().equals(otherPair.second(), true))
                    return false;
            }

            return true;
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

                final SvLinkedPair otherPair = other.getLinkedPairs().get(j);

                boolean linksMatch = pair.first().equals(otherPair.first(), true)
                    && pair.second().equals(otherPair.second(), true);

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
            else if(pair.locationType() == LOCATION_TYPE_REMOTE)
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
