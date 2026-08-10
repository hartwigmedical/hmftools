package com.hartwig.hmftools.bamtools.synthetic;

import static java.lang.Math.max;
import static java.lang.Math.min;
import static java.lang.String.format;

import java.util.List;
import java.util.stream.Collectors;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.region.ChrBaseRegion;

public class RegionData extends ChrBaseRegion
{
    private final List<RegionType> mTypes;
    private double mMaxCopyNumber;

    public RegionData(final String chromosome, final int posStart, final int posEnd, final RegionType regionType)
    {
        super(chromosome, posStart, posEnd);
        mTypes = Lists.newArrayList(regionType);
        mMaxCopyNumber = 0;
    }

    public List<RegionType> types() { return mTypes; }

    public void addRegionTypes(final List<RegionType> types)
    {
        types.stream().filter(x -> !mTypes.contains(x)).forEach(x -> mTypes.add(x));
    }

    public double copyNumber() { return mMaxCopyNumber; }
    public void setCopyNumber(double copyNumber) { mMaxCopyNumber = max(mMaxCopyNumber, copyNumber); }

    public void merge(final RegionData other)
    {
        setStart(min(start(), other.start()));
        setEnd(max(end(), other.end()));
        addRegionTypes(other.types());
        setCopyNumber(other.copyNumber());
    }

    public String toString()
    {
        return format("%s: %s", super.toString(), mTypes.stream().map(x -> x.toString()).collect(Collectors.joining(";")));
    }
}
