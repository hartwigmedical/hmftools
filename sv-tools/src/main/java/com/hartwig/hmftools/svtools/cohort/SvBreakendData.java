package com.hartwig.hmftools.svtools.cohort;

import com.hartwig.hmftools.common.sv.StructuralVariantType;

public class SvBreakendData
{
    public final String SampleId;
    public final int SvId;
    public final int Position;
    public final boolean IsStart;
    public final byte Orientation;
    public final StructuralVariantType Type;

    public SvBreakendData(final String sampleId, final int svId, final int position, final boolean isStart, final byte orientation,
            final StructuralVariantType type)
    {
        SampleId = sampleId;
        SvId = svId;
        Position = position;
        IsStart = isStart;
        Orientation = orientation;
        Type = type;
    }
}
