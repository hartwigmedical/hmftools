package com.hartwig.hmftools.common.stats;

import org.jetbrains.annotations.NotNull;

public class PValueResult implements Comparable<PValueResult>
{
    public final String Id;
    public final double PValue;
    public int Rank;
    public double QValue;

    public PValueResult(final String id, final double pValue)
    {
        Id = id;
        PValue = pValue;
        Rank = 0;
        QValue = 0;
    }

    public int compareTo(@NotNull PValueResult other)
    {
        if(PValue == other.PValue)
            return 0;

        return PValue < other.PValue ? -1 : 1;
    }

    public String toString() { return String.format("%s: p(%4.3e) q(%4.3e) rank(%d)", Id, PValue, QValue, Rank); }
}
