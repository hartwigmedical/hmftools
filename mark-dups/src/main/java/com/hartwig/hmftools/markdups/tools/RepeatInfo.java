package com.hartwig.hmftools.markdups.tools;

import java.util.List;
import java.util.stream.Collectors;

import com.hartwig.hmftools.common.region.BaseRegion;
import com.hartwig.hmftools.common.region.ChrBaseRegion;

public class RepeatInfo
{
    private final List<Character> mRepeatBases;
    private final ChrBaseRegion mRepeatRegion;
    private final List<BaseRegion> mContainedWithinBedRegions;
    private final List<BaseRegion> mPartialOverlapsWithBedRegions;
    private boolean mBedAnnotations;

    public RepeatInfo(final List<Character> repeatBases, final ChrBaseRegion repeatRegion, final List<BaseRegion> containedWithinBedRegions,
            final List<BaseRegion> partialOverlapsWithBedRegions)
    {
        mRepeatBases = repeatBases;
        mRepeatRegion = repeatRegion;
        mContainedWithinBedRegions = containedWithinBedRegions;
        mPartialOverlapsWithBedRegions = partialOverlapsWithBedRegions;
        mBedAnnotations = true;
    }

    public RepeatInfo(final List<Character> repeatBases, final ChrBaseRegion repeatRegion)
    {
        this(repeatBases, repeatRegion, null, null);
        mBedAnnotations = false;
    }

    public String getTSVFragment()
    {
        String repeatBasesStr = String.valueOf(mRepeatBases.get(0));
        if(!mRepeatBases.get(0).equals(mRepeatBases.get(1)))
        {
            repeatBasesStr = repeatBasesStr + mRepeatBases.get(1);
        }

        String output =
                String.format("%s\t%d\t%d\t%s", mRepeatRegion.Chromosome, mRepeatRegion.start(), mRepeatRegion.end(), repeatBasesStr);
        if(!mBedAnnotations)
        {
            return output;
        }

        final String containedWithinRegionsStr =
                mContainedWithinBedRegions.stream().map(x -> x.toString()).collect(Collectors.joining(","));
        final String partialOverlapsWithRegionsStr =
                mPartialOverlapsWithBedRegions.stream().map(x -> x.toString()).collect(Collectors.joining(","));
        return output + '\t' + containedWithinRegionsStr + '\t' + partialOverlapsWithRegionsStr;
    }
}
