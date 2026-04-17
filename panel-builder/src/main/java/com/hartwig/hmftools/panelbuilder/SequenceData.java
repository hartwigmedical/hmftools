package com.hartwig.hmftools.panelbuilder;

public record SequenceData(
        byte[] bases,
        boolean isNormal,
        int gcCount
)
{
    public String baseString()
    {
        return new String(bases);
    }

    public double gcContent()
    {
        return gcCount / (double) bases.length;
    }
}
