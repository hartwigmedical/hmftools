package com.hartwig.hmftools.common.aligner;

import com.google.common.base.Preconditions;

class AlignerTraits
{
    // scores
    protected final int mMatchScore;
    protected final int mMismatchScore;
    protected final int mGapOpeningScore;
    protected final int mGapExtensionScore;

    protected boolean mLogWorkMatrix = false;

    public AlignerTraits(int matchScore, int mismatchScore, int gapOpeningScore, int gapExtensionScore)
    {
        Preconditions.checkArgument(matchScore > 0);
        Preconditions.checkArgument(mismatchScore <= 0);
        Preconditions.checkArgument(gapOpeningScore <= 0);
        Preconditions.checkArgument(gapExtensionScore <= 0);
        this.mMatchScore = matchScore;
        this.mMismatchScore = mismatchScore;
        this.mGapOpeningScore = gapOpeningScore;
        this.mGapExtensionScore = gapExtensionScore;
    }

    public void setLogWorkMatrix(boolean b)
    {
        mLogWorkMatrix = b;
    }
}
