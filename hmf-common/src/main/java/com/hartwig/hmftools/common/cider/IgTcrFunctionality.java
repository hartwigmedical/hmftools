package com.hartwig.hmftools.common.cider;

// see https://www.imgt.org/IMGTScientificChart/SequenceDescription/IMGTfunctionality.html
public enum IgTcrFunctionality
{
    FUNCTIONAL,
    ORF, // open reading frame
    PSEUDOGENE; // pseudogene

    public String toCode()
    {
        return switch(this)
        {
            case FUNCTIONAL -> "F";
            case ORF -> "ORF";
            case PSEUDOGENE -> "P";
        };
    }

    public static IgTcrFunctionality fromCode(String code)
    {
        return switch(code)
        {
            case "F" -> FUNCTIONAL;
            case "ORF" -> ORF;
            case "P" -> PSEUDOGENE;
            default -> throw new IllegalArgumentException("invalid IgTcrFunctionality code: " + code);
        };
    }
}
