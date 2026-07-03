package com.hartwig.hmftools.tars.liftback.overhang;

// Counters for the overhang gate; aggregated across workers and logged at end-of-run. Junctions collapsed
// (leading/trailing) come from the peel; reclaimed records/bases from the standalone softclip reclaim; alts
// dropped from XA-alt placements the peel collapsed to a contiguous alignment.
public class OverhangGateStatistics
{
    private long mCollapsedLeading;
    private long mCollapsedTrailing;
    private long mRecordsReclaimed;
    private long mBasesReclaimedLead;
    private long mBasesReclaimedTrail;
    private long mAltsDropped;

    public void countCollapsedLeading() { ++mCollapsedLeading; }

    public void countCollapsedTrailing() { ++mCollapsedTrailing; }

    public void countReclaimed(final long leadBases, final long trailBases)
    {
        ++mRecordsReclaimed;
        mBasesReclaimedLead += leadBases;
        mBasesReclaimedTrail += trailBases;
    }

    public void countAltDropped() { ++mAltsDropped; }

    public long collapsedLeading() { return mCollapsedLeading; }

    public long collapsedTrailing() { return mCollapsedTrailing; }

    public long recordsReclaimed() { return mRecordsReclaimed; }

    public long basesReclaimedLead() { return mBasesReclaimedLead; }

    public long basesReclaimedTrail() { return mBasesReclaimedTrail; }

    public long altsDropped() { return mAltsDropped; }

    // Snapshot / restore so a discarded provisional mate decision can be rolled back without double-counting
    // (see LiftBackGroupProcessor.processNameGroup).
    public long[] snapshot()
    {
        return new long[] {
                mCollapsedLeading, mCollapsedTrailing, mRecordsReclaimed,
                mBasesReclaimedLead, mBasesReclaimedTrail, mAltsDropped };
    }

    public void restore(final long[] snapshot)
    {
        mCollapsedLeading = snapshot[0];
        mCollapsedTrailing = snapshot[1];
        mRecordsReclaimed = snapshot[2];
        mBasesReclaimedLead = snapshot[3];
        mBasesReclaimedTrail = snapshot[4];
        mAltsDropped = snapshot[5];
    }
}
