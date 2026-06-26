package com.hartwig.hmftools.tars.liftback.tailextend;

// Counters for the tail-extension pass; logged at end-of-run.
public class TailExtensionStatistics
{
    private int mRecordsEvaluated;
    private int mRecordsExtended;
    private int mBasesExtendedLead;
    private int mBasesExtendedTrail;
    private int mSkippedNoRef;
    private int mSkippedForJunctionGuard;
    private int mSkippedComplexShape;
    private int mRejectedTooManyMismatches;

    public void countEvaluated()
    {
        ++mRecordsEvaluated;
    }

    public void countExtended(final int leadBases, final int trailBases)
    {
        ++mRecordsExtended;
        mBasesExtendedLead += leadBases;
        mBasesExtendedTrail += trailBases;
    }

    public void countSkippedNoRef()
    {
        ++mSkippedNoRef;
    }

    public void countSkippedForJunctionGuard()
    {
        ++mSkippedForJunctionGuard;
    }

    public void countSkippedComplexShape()
    {
        ++mSkippedComplexShape;
    }

    public void countRejectedTooManyMismatches()
    {
        ++mRejectedTooManyMismatches;
    }

    public int recordsEvaluated() { return mRecordsEvaluated; }

    public int recordsExtended() { return mRecordsExtended; }

    public int basesExtendedLead() { return mBasesExtendedLead; }

    public int basesExtendedTrail() { return mBasesExtendedTrail; }

    public int skippedNoRef() { return mSkippedNoRef; }

    public int skippedForJunctionGuard() { return mSkippedForJunctionGuard; }

    public int skippedComplexShape() { return mSkippedComplexShape; }

    public int rejectedTooManyMismatches() { return mRejectedTooManyMismatches; }

    // Snapshot / restore the counters so a discarded provisional mate decision can be rolled back without
    // double-counting (see LiftBackGroupProcessor.processNameGroup).
    public int[] snapshot()
    {
        return new int[] {
                mRecordsEvaluated, mRecordsExtended, mBasesExtendedLead, mBasesExtendedTrail,
                mSkippedNoRef, mSkippedForJunctionGuard, mSkippedComplexShape, mRejectedTooManyMismatches };
    }

    public void restore(final int[] snapshot)
    {
        mRecordsEvaluated = snapshot[0];
        mRecordsExtended = snapshot[1];
        mBasesExtendedLead = snapshot[2];
        mBasesExtendedTrail = snapshot[3];
        mSkippedNoRef = snapshot[4];
        mSkippedForJunctionGuard = snapshot[5];
        mSkippedComplexShape = snapshot[6];
        mRejectedTooManyMismatches = snapshot[7];
    }
}
