package com.hartwig.hmftools.purple.tools;

import static com.hartwig.hmftools.common.region.BaseRegion.positionsOverlap;

import java.util.List;

import com.google.common.base.Preconditions;
import com.hartwig.hmftools.common.gene.ExonData;

public class GeneAmplification
{
    private final int mTotalExonCount;
    private final int mOverlappingExonCount;
    private final boolean mOverlapsFirstExon;
    private final boolean mOverlapsLastExon;

    GeneAmplification(List<ExonData> allExons, int start, int end)
    {
        Preconditions.checkArgument(!allExons.isEmpty());
        Preconditions.checkArgument(start < end);

        mTotalExonCount = allExons.size();
        List<ExonData> overlappingExons = allExons.stream().filter(x -> positionsOverlap(x.Start, x.End, start, end)).toList();
        mOverlappingExonCount = overlappingExons.size();

        if(mOverlappingExonCount > 0)
        {
            mOverlapsFirstExon = overlappingExons.get(0).equals(allExons.get(0));
            mOverlapsLastExon = overlappingExons.get(mOverlappingExonCount - 1).equals(allExons.get(mTotalExonCount - 1));
        }
        else
        {
            mOverlapsFirstExon = false;
            mOverlapsLastExon = false;
        }
    }

    public boolean isCompleteAmplification()
    {
        return mOverlappingExonCount == mTotalExonCount;
    }

    public boolean isOfInterest()
    {
        if(mOverlappingExonCount == 0)
        {
            return false;
        }
        if(isCompleteAmplification())
        {
            return true;
        }
        return !mOverlapsFirstExon && !mOverlapsLastExon;
    }

    public boolean isTailAmplification()
    {
        return mOverlappingExonCount > 0 && !isCompleteAmplification() && mOverlapsLastExon;
    }

    public boolean isHeadAmplification()
    {
        return mOverlappingExonCount > 0 && !isCompleteAmplification() && mOverlapsFirstExon;
    }

    public int numberOfAffectedExons()
    {
        return mOverlappingExonCount;
    }
}
