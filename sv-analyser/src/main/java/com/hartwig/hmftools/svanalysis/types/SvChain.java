package com.hartwig.hmftools.svanalysis.types;

import java.util.List;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.svanalysis.analysis.SvCluster;

public class SvChain {

    private int mChainId;
    private List<SvClusterData> mSVs;
    private List<SvLinkedPair> mLinkedPairs;

    private List<Integer> mTILengths;

    public SvChain(int chainId)
    {
        mChainId = chainId;
        mSVs = Lists.newArrayList();
        mLinkedPairs = Lists.newArrayList();
        mTILengths = Lists.newArrayList();
    }

    public int getId() { return mChainId; }

    public int getCount() { return mSVs.size(); }
    // public List<SvClusterData> getSVs() { return mSVs; }

    public List<SvLinkedPair> getLinkedPairs() { return mLinkedPairs; }

    public int getLinkCount() { return mLinkedPairs.size(); }

    public void addVariant(final SvClusterData variant, int tiLength)
    {
        mSVs.add(variant);

        if(tiLength > 0)
        {
            mTILengths.add(tiLength);
        }
    }

    public void addLinkedPair(final SvLinkedPair pair, boolean addToStart)
    {
        if(mLinkedPairs.isEmpty() || !addToStart)
            mLinkedPairs.add(pair);
        else
            mLinkedPairs.add(0, pair);
    }

    public SvClusterData getFirstSV() { return mLinkedPairs.get(0).first(); }
    public boolean firstUnlinkedOnStart() { return !mLinkedPairs.get(0).firstLinkOnStart(); }

    public SvClusterData getLastSV() { return mLinkedPairs.get(mLinkedPairs.size() - 1).second(); }
    public boolean lastUnlinkedOnStart() { return !mLinkedPairs.get(0).secondLinkOnStart(); }

    public int getTICount()
    {
        int tiCount = 0;
        for(final SvLinkedPair pair : mLinkedPairs)
        {
            if(pair.linkType() == SvLinkedPair.LINK_TYPE_TI)
                ++tiCount;
        }

        return tiCount;
    }

    public int getDBCount()
    {
        int dbCount = 0;
        for(final SvLinkedPair pair : mLinkedPairs)
        {
            if(pair.linkType() == SvLinkedPair.LINK_TYPE_DB)
                ++dbCount;
        }

        return dbCount;
    }

    public boolean hasOnlyShortTIs(long maxLength)
    {
        for(final SvLinkedPair pair : mLinkedPairs)
        {
            if(pair.linkType() != SvLinkedPair.LINK_TYPE_TI || pair.length() > maxLength )
                return false;
        }

        return true;
    }

    public int getSvIndex(final SvClusterData var)
    {
        for(int index = 0; index < mLinkedPairs.size(); ++index)
        {
            final SvLinkedPair pair = mLinkedPairs.get(index);

            if(pair.first().equals(var) || pair.second().equals(var))
            {
                return index;
            }
        }

        return -1;
    }

    public int getSvIndex(final SvClusterData var, boolean matchStart)
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

    public boolean hasLinkedPair(final SvLinkedPair pair)
    {
        return mLinkedPairs.contains(pair);
    }

    public boolean hasVariantBE(final SvClusterData var, boolean useStart)
    {
        for(final SvLinkedPair pair : mLinkedPairs)
        {
            if(pair.hasVariantBE(var, useStart))
                return true;
        }

        return false;
    }


}
