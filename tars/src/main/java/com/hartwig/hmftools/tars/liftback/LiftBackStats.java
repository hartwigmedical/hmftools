package com.hartwig.hmftools.tars.liftback;

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
    public long SupplementaryMerges;
    public long SupplementariesAbsorbed;

    public void add(final LiftBackStats other)
    {
        RecordsSeen += other.RecordsSeen;
        PrimariesSeen += other.PrimariesSeen;
        LiftFailed += other.LiftFailed;
        UnmappedExcludedRegion += other.UnmappedExcludedRegion;
        UnmappedOverCap += other.UnmappedOverCap;
        UnmappedLowAlignmentScore += other.UnmappedLowAlignmentScore;
        MergeableSupplementaries += other.MergeableSupplementaries;
        SupplementaryMerges += other.SupplementaryMerges;
        SupplementariesAbsorbed += other.SupplementariesAbsorbed;
    }
}
