package com.hartwig.hmftools.esvee.models;

public class AssemblyClassification
{
    public final AssemblyClassificationType Type;

    /** The length of this event, if applicable. 0 otherwise */
    public final int Length;

    public AssemblyClassification(final AssemblyClassificationType type, final int length)
    {
        Type = type;
        Length = length;
    }

    @Override
    public String toString()
    {
        final boolean showLength = Type != AssemblyClassificationType.TRANSLOCATION
                && Type != AssemblyClassificationType.INVERSION
                && Type != AssemblyClassificationType.UNKNOWN;
        return Type.name() + (showLength ? " " + Length : "");
    }
}
