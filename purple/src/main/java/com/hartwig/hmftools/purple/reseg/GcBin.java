package com.hartwig.hmftools.purple.reseg;

import static com.hartwig.hmftools.common.utils.Doubles.round;
import static com.hartwig.hmftools.purple.PurpleConstants.RESEG_GC_BIN_MAX;
import static com.hartwig.hmftools.purple.PurpleConstants.RESEG_GC_BIN_MIN;
import static com.hartwig.hmftools.purple.PurpleConstants.RESEG_GC_BIN_WIDTH;

import java.util.ArrayList;
import java.util.List;

public final class GcBin
{
    public final double Lower;
    public final double Upper;

    public GcBin(final double lower, final double upper)
    {
        Lower = lower;
        Upper = upper;
    }

    public boolean contains(double gcContent)
    {
        return gcContent > Lower && gcContent <= Upper;
    }

    public String toString()
    {
        return String.format("(%.2f-%.2f]", Lower, Upper);
    }

    public static List<GcBin> createdBins()
    {
        List<GcBin> bins = new ArrayList<>();
        int binCount = (int)Math.round((RESEG_GC_BIN_MAX - RESEG_GC_BIN_MIN) / RESEG_GC_BIN_WIDTH);

        for(int i = 0; i < binCount; i++)
        {
            double lower = round(RESEG_GC_BIN_MIN + i * RESEG_GC_BIN_WIDTH, 2);
            double upper = round(lower + RESEG_GC_BIN_WIDTH, 2);
            bins.add(new GcBin(lower, upper));
        }

        return bins;
    }
}
