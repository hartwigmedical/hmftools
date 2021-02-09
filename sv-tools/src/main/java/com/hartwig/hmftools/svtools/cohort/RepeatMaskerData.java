package com.hartwig.hmftools.svtools.cohort;

import static java.lang.Math.abs;

import com.hartwig.hmftools.common.utils.sv.BaseRegion;

public class RepeatMaskerData
{
    public final int RmId;
    public final BaseRegion Region;
    public final byte Strand;

    public RepeatMaskerData(final int rmId, final BaseRegion region, final byte strand)
    {
        RmId = rmId;
        Region = region;
        Strand = strand;
    }

    public String toString() { return String.format("%d: region(%s) strand(%d)", RmId, Region, Strand); }

    public boolean isProximateElement(final RepeatMaskerData other, final int maxDistance)
    {
        if(abs(Region.start() - other.Region.end()) <= maxDistance)
            return true;

        if(abs(other.Region.start() - Region.end()) <= maxDistance)
            return true;

        return false;
    }

}
