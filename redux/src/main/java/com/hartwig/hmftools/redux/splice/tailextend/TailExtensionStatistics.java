package com.hartwig.hmftools.redux.splice.tailextend;

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
}
