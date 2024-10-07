package com.hartwig.hmftools.chord.indel;

public enum ContextType
{
    REPEAT(".rep.len."),
    MICROHOMOLOGY(".mh.bimh."),
    NONE(".none.len.");

    public final String Suffix;

    ContextType(String suffix)
    {
        Suffix = suffix;
    }
}
