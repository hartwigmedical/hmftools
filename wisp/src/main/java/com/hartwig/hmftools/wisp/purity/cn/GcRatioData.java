package com.hartwig.hmftools.wisp.purity.cn;

public class GcRatioData implements Comparable<GcRatioData>
{
    public final int Position;
    public final double TumorGcRatio;

    public GcRatioData(final int position, final double tumorGcRatio)
    {
        Position = position;
        TumorGcRatio = tumorGcRatio;
    }

    @Override
    public int compareTo(final GcRatioData other)
    {
        if(TumorGcRatio == other.TumorGcRatio)
            return 0;

        return TumorGcRatio < other.TumorGcRatio ? -1 : 1;
    }
}
