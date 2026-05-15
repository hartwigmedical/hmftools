package com.hartwig.hmftools.redux.splice;

import java.util.HashMap;
import java.util.Map;

// in-memory cache of lifted primary alignment info, keyed by read name. Built in pass 1 of SpliceLiftBack
// (one entry per paired primary record) and consulted in pass 2 to patch each emitted record's mate fields
// with its partner's lifted coordinates.
//
// Holds only primary alignments. Supplementary and secondary records do not populate mate fields in other
// records, so they are not cached.
public class LiftedMateInfoCache
{
    private final Map<String, ReadPairLiftedMateInfo> mLiftedMateInfoByReadName;

    public LiftedMateInfoCache()
    {
        mLiftedMateInfoByReadName = new HashMap<>();
    }

    public void recordPrimaryAlignment(final String readName, final boolean firstOfPair, final LiftedMateInfo info)
    {
        ReadPairLiftedMateInfo pairInfo = mLiftedMateInfoByReadName.computeIfAbsent(readName, k -> new ReadPairLiftedMateInfo());
        if(firstOfPair)
            pairInfo.FirstInPair = info;
        else
            pairInfo.SecondInPair = info;
    }

    // returns the lifted info for the partner of the given read (i.e. R2's info when looking up R1, and vice versa),
    // or null if the partner's primary was not seen in pass 1.
    public LiftedMateInfo getPartnerMateInfo(final String readName, final boolean firstOfPair)
    {
        ReadPairLiftedMateInfo pairInfo = mLiftedMateInfoByReadName.get(readName);
        if(pairInfo == null)
            return null;
        return firstOfPair ? pairInfo.SecondInPair : pairInfo.FirstInPair;
    }

    // returns the lifted info for this read's own primary (i.e. R1's primary when given a non-primary R1),
    // used to mirror primary coords onto a supplementary record whose own lift failed.
    public LiftedMateInfo getOwnPrimary(final String readName, final boolean firstOfPair)
    {
        ReadPairLiftedMateInfo pairInfo = mLiftedMateInfoByReadName.get(readName);
        if(pairInfo == null)
            return null;
        return firstOfPair ? pairInfo.FirstInPair : pairInfo.SecondInPair;
    }

    public int size()
    {
        return mLiftedMateInfoByReadName.size();
    }

    private static class ReadPairLiftedMateInfo
    {
        public LiftedMateInfo FirstInPair;
        public LiftedMateInfo SecondInPair;
    }
}
