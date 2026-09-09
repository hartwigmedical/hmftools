package com.hartwig.hmftools.tars.liftback;

import com.hartwig.hmftools.tars.liftback.features.SupplementaryMerger.RejectReason;

// Counters for one worker, summed across workers for the run summary. The three deliberate unmapping routes are normal
// outcomes, kept apart from LiftFailed which alone signals a sidecar/FASTA mismatch.
public class LiftBackStats
{
    public long RecordsSeen;
    public long PrimariesSeen;
    public long LiftFailed;
    public long UnmappedExcludedRegion;
    public long UnmappedOverCap;
    public long UnmappedLowAlignmentScore;
    public long MergeableSupplementaries;
    public long SupplementaryMergeAttempts;
    public long SupplementaryMergeSuccessfulCandidates;
    public long SupplementaryMerges;
    public long SupplementariesAbsorbed;
    private final long[] mSupplementaryMergeRejections = new long[RejectReason.values().length];

    public void recordSupplementaryMergeRejection(final RejectReason reason)
    {
        ++mSupplementaryMergeRejections[reason.ordinal()];
    }

    public long supplementaryMergeRejections(final RejectReason reason)
    {
        return mSupplementaryMergeRejections[reason.ordinal()];
    }

    public void add(final LiftBackStats other)
    {
        RecordsSeen += other.RecordsSeen;
        PrimariesSeen += other.PrimariesSeen;
        LiftFailed += other.LiftFailed;
        UnmappedExcludedRegion += other.UnmappedExcludedRegion;
        UnmappedOverCap += other.UnmappedOverCap;
        UnmappedLowAlignmentScore += other.UnmappedLowAlignmentScore;
        MergeableSupplementaries += other.MergeableSupplementaries;
        SupplementaryMergeAttempts += other.SupplementaryMergeAttempts;
        SupplementaryMergeSuccessfulCandidates += other.SupplementaryMergeSuccessfulCandidates;
        SupplementaryMerges += other.SupplementaryMerges;
        SupplementariesAbsorbed += other.SupplementariesAbsorbed;
        for(RejectReason reason : RejectReason.values())
        {
            mSupplementaryMergeRejections[reason.ordinal()] += other.supplementaryMergeRejections(reason);
        }
    }
}
