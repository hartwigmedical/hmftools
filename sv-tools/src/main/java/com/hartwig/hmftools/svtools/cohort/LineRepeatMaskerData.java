package com.hartwig.hmftools.svtools.cohort;

import com.hartwig.hmftools.common.utils.sv.SvRegion;

public class LineRepeatMaskerData
{
    public final int RmId;
    public final SvRegion Region;
    public final byte Strand;

    public LineRepeatMaskerData(final int rmId, final SvRegion region, final byte strand)
    {
        RmId = rmId;
        Region = region;
        Strand = strand;
    }
}
