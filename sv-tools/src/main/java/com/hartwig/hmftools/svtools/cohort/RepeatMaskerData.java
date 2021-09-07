package com.hartwig.hmftools.svtools.cohort;

import static java.lang.Math.abs;

import com.hartwig.hmftools.common.utils.sv.ChrBaseRegion;

public class RepeatMaskerData
{
    public final int RmId;
    public final ChrBaseRegion Region;
    public final byte Strand;
    public final String ClassFamily;
    public final String Repeat;

    public RepeatMaskerData(final int rmId, final ChrBaseRegion region, final byte strand)
    {
        RmId = rmId;
        Region = region;
        Strand = strand;
        ClassFamily = "";
        Repeat = "";
    }

    public RepeatMaskerData(final int rmId, final ChrBaseRegion region, final byte strand, final String classFamily, final String repeat)
    {
        RmId = rmId;
        Region = region;
        Strand = strand;
        ClassFamily = classFamily;
        Repeat = repeat;
    }

    public String toString() { return String.format("%d: region(%s) strand(%d) repeat(%s)",
            RmId, Region, Strand, Repeat); }

    public boolean isProximateElement(final RepeatMaskerData other, final int maxDistance)
    {
        if(abs(Region.start() - other.Region.end()) <= maxDistance)
            return true;

        if(abs(other.Region.start() - Region.end()) <= maxDistance)
            return true;

        return false;
    }

    public String subRepeat(int length) { return Repeat.length() < length ? Repeat : Repeat.substring(0, length); }

}
