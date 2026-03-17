package com.hartwig.hmftools.common.amber;

public enum AmberBase
{
    G,
    A,
    T,
    C,
    N;

    public AmberBase complement()
    {
        return switch(this)
        {
            case A -> T;
            case T -> A;
            case C -> G;
            case G -> C;
            default -> N;
        };
    }
}
