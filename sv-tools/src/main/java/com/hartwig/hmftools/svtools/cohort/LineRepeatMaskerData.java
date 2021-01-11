package com.hartwig.hmftools.svtools.cohort;

import com.hartwig.hmftools.common.utils.sv.BaseRegion;

public class LineRepeatMaskerData
{
    public final int RmId;
    public final BaseRegion Region;
    public final byte Strand;

    public LineRepeatMaskerData(final int rmId, final BaseRegion region, final byte strand)
    {
        RmId = rmId;
        Region = region;
        Strand = strand;
    }
}
