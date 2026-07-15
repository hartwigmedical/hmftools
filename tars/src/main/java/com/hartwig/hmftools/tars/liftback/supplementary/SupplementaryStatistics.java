package com.hartwig.hmftools.tars.liftback.supplementary;

import java.util.EnumMap;
import java.util.Map;

// Per-decision counters for the supplementary resolver. Aggregated across the run and emitted into the
// liftback summary so we can see which gate is rejecting most candidates and how many chain
// merges happened.
public class SupplementaryStatistics
{
    private int mCandidatesEvaluated;
    private int mMergedTotal;
    private int mSuppClampApplied;
    private final int[] mChainDepthCounts;  // index = depth (1 = single merge, 2 = chain of two, ...)
    private final Map<SupplementaryRejectReason, Integer> mRejections;

    public SupplementaryStatistics(final int maxChainDepth)
    {
        mChainDepthCounts = new int[Math.max(maxChainDepth, 1) + 1];
        mRejections = new EnumMap<>(SupplementaryRejectReason.class);
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

    public void countReject(final SupplementaryRejectReason reason)
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

    public int rejectCount(final SupplementaryRejectReason reason)
    {
        return mRejections.getOrDefault(reason, 0);
    }

    // Snapshot / restore all counters so a discarded provisional mate decision can be rolled back without
    // double-counting (see LiftBackGroupProcessor.processNameGroup).
    public int[] snapshot()
    {
        SupplementaryRejectReason[] reasons = SupplementaryRejectReason.values();
        int[] snapshot = new int[3 + mChainDepthCounts.length + reasons.length];
        snapshot[0] = mCandidatesEvaluated;
        snapshot[1] = mMergedTotal;
        snapshot[2] = mSuppClampApplied;
        System.arraycopy(mChainDepthCounts, 0, snapshot, 3, mChainDepthCounts.length);
        for(int i = 0; i < reasons.length; ++i)
        {
            snapshot[3 + mChainDepthCounts.length + i] = mRejections.getOrDefault(reasons[i], 0);
        }
        return snapshot;
    }

    public void restore(final int[] snapshot)
    {
        mCandidatesEvaluated = snapshot[0];
        mMergedTotal = snapshot[1];
        mSuppClampApplied = snapshot[2];
        System.arraycopy(snapshot, 3, mChainDepthCounts, 0, mChainDepthCounts.length);
        SupplementaryRejectReason[] reasons = SupplementaryRejectReason.values();
        mRejections.clear();
        for(int i = 0; i < reasons.length; ++i)
        {
            int count = snapshot[3 + mChainDepthCounts.length + i];
            if(count != 0)
            {
                mRejections.put(reasons[i], count);
            }
        }
    }
}
