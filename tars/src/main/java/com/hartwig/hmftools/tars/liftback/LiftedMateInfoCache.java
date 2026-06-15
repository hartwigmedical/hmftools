package com.hartwig.hmftools.tars.liftback;

import java.util.HashMap;
import java.util.Map;

// Cache of lifted primary alignment info keyed by read name. Built in pass 1, consumed in pass 2 to patch
// mate fields (RNEXT/PNEXT/TLEN). Supplementary/secondary records are not cached — they don't populate mates.
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

    public LiftedMateInfo getPartnerMateInfo(final String readName, final boolean firstOfPair)
    {
        ReadPairLiftedMateInfo pairInfo = mLiftedMateInfoByReadName.get(readName);
        if(pairInfo == null)
            return null;
        return firstOfPair ? pairInfo.SecondInPair : pairInfo.FirstInPair;
    }

    // used to mirror primary coords onto a supplementary whose own lift failed
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
