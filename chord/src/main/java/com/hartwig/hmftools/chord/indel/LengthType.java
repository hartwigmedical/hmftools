package com.hartwig.hmftools.chord.indel;

public enum LengthType
{
    INDEL_LENGTH("len"), // Same as repeat unit length
    HOM_BASES_COUNT("bimh"); // Bases in microhomology

    public final String Abbreviation;

    LengthType(String abbreviation)
    {
        Abbreviation = abbreviation;
    }

    public LengthType from(ContextType contextType)
    {
        return (contextType == ContextType.MICROHOMOLOGY) ? HOM_BASES_COUNT : INDEL_LENGTH;
    }
}
