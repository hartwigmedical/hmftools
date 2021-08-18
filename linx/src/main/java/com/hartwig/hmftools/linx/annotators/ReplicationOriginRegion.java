package com.hartwig.hmftools.linx.annotators;

import com.hartwig.hmftools.common.utils.sv.ChrBaseRegion;

public class ReplicationOriginRegion
{
    public final ChrBaseRegion Region;
    public final double OriginValue;

    public ReplicationOriginRegion(final String chromosome, int start, int end, double originValue)
    {
        Region = new ChrBaseRegion(chromosome, start, end);
        OriginValue = originValue;
    }
}
