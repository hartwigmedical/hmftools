package com.hartwig.hmftools.esvee.assembly.types;

public enum BaseType
{
    A('A'),
    C('C'),
    G('G'),
    T('T'),
    OTHER('X');

    public final byte Byte;
    public final char Char;

    BaseType(final char asChar)
    {
        Char = asChar;
        Byte = (byte)asChar;
    }

    public static BaseType from(final byte base)
    {
        if(base == A.Byte) return A;
        if(base == C.Byte) return C;
        if(base == G.Byte) return G;
        if(base == T.Byte) return T;
        return OTHER;
    }
}
