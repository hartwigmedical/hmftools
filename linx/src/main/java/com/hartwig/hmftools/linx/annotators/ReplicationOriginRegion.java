package com.hartwig.hmftools.linx.annotators;

import com.hartwig.hmftools.common.utils.sv.BaseRegion;

public class ReplicationOriginRegion
{
    public final BaseRegion Region;
    public final double OriginValue;

    public ReplicationOriginRegion(final String chromosome, int start, int end, double originValue)
    {
        Region = new BaseRegion(chromosome, start, end);
        OriginValue = originValue;
    }
}
