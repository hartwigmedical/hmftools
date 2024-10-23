package com.hartwig.hmftools.esvee.alignment;

import static java.lang.String.format;

public class MdTagElement
{
    public final MdTagType Type;
    public final int Length;
    public final char Base;

    public MdTagElement(final MdTagType type, final int length, final char base)
    {
        Type = type;
        Length = length;
        Base = base;
    }

    // convenience
    public boolean isMatch() { return Type == MdTagType.MATCH; }
    public boolean isIndel() { return Type == MdTagType.DEL; }
    public boolean isSnv() { return Type == MdTagType.SNV; }

    public String toString()
    {
        if(Type == MdTagType.SNV)
            return String.valueOf(Base);
        else if(Type == MdTagType.MATCH)
            return format("%dM", Length);
        else
            return format("%dD", Length);
    }
}
