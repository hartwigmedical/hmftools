package com.hartwig.hmftools.redux.merge.repeatfinder;

import java.util.List;
import java.util.stream.Collectors;

import com.hartwig.hmftools.common.region.BaseRegion;
import com.hartwig.hmftools.common.region.ChrBaseRegion;

public class RepeatInfo
{
    private final ChrBaseRegion mRepeatRegion;
    private final String mRepeatBasesStr;
    private final List<BaseRegion> mContainedWithinBedRegions;
    private final List<BaseRegion> mPartialOverlapsWithBedRegions;
    private boolean mBedAnnotations;

    public RepeatInfo(final ChrBaseRegion repeatRegion, final String repeatBasesStr, final List<BaseRegion> containedWithinBedRegions,
            final List<BaseRegion> partialOverlapsWithBedRegions)
    {
        mRepeatRegion = repeatRegion;
        mRepeatBasesStr = repeatBasesStr;
        mContainedWithinBedRegions = containedWithinBedRegions;
        mPartialOverlapsWithBedRegions = partialOverlapsWithBedRegions;
        mBedAnnotations = true;
    }

    public RepeatInfo(final ChrBaseRegion repeatRegion, final String repeatBasesStr)
    {
        this(repeatRegion, repeatBasesStr, null, null);
        mBedAnnotations = false;
    }

    public String getTSVFragment()
    {
        String output =
                String.format("%s\t%d\t%d\t%s", mRepeatRegion.Chromosome, mRepeatRegion.start(), mRepeatRegion.end(), mRepeatBasesStr);
        if(!mBedAnnotations)
        {
            return output;
        }

        String containedWithinRegionsStr =
                mContainedWithinBedRegions.stream().map(x -> x.toString()).collect(Collectors.joining(","));
        String partialOverlapsWithRegionsStr =
                mPartialOverlapsWithBedRegions.stream().map(x -> x.toString()).collect(Collectors.joining(","));
        return output + '\t' + containedWithinRegionsStr + '\t' + partialOverlapsWithRegionsStr;
    }
}
