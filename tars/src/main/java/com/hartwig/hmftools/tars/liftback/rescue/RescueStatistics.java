package com.hartwig.hmftools.tars.liftback.rescue;

import java.util.EnumMap;
import java.util.Map;

// Per-decision counters for the rescue resolver. Aggregated across the run and emitted into the
// liftback summary so we can see which gate is rejecting most candidates and how many chain
// merges happened.
public class RescueStatistics
{
    private int mCandidatesEvaluated;
    private int mMergedTotal;
    private int mSuppClampApplied;
    private final int[] mChainDepthCounts;  // index = depth (1 = single merge, 2 = chain of two, ...)
    private final Map<RescueRejectReason, Integer> mRejections;

    public RescueStatistics(final int maxChainDepth)
    {
        mChainDepthCounts = new int[Math.max(maxChainDepth, 1) + 1];
        mRejections = new EnumMap<>(RescueRejectReason.class);
    }

    public void countCandidate()
    {
        ++mCandidatesEvaluated;
    }

    public void countMergedChain(final int chainDepth)
    {
        ++mMergedTotal;
        int slot = Math.min(chainDepth, mChainDepthCounts.length - 1);
        ++mChainDepthCounts[slot];
    }

    public void countReject(final RescueRejectReason reason)
    {
        mRejections.merge(reason, 1, Integer::sum);
    }

    public void countSuppClamp()
    {
        ++mSuppClampApplied;
    }

    public int suppClampApplied()
    {
        return mSuppClampApplied;
    }

    public int candidatesEvaluated()
    {
        return mCandidatesEvaluated;
    }

    public int mergedTotal()
    {
        return mMergedTotal;
    }

    public int mergedAtChainDepth(final int depth)
    {
        if(depth <= 0 || depth >= mChainDepthCounts.length)
        {
            return 0;
        }
        return mChainDepthCounts[depth];
    }

    public int rejectCount(final RescueRejectReason reason)
    {
        return mRejections.getOrDefault(reason, 0);
    }
}
