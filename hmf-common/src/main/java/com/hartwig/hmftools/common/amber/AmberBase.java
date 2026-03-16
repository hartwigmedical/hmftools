package com.hartwig.hmftools.common.amber;

public enum AmberBase
{
    G
            {
                @Override
                public AmberBase complement()
                {
                    return C;
                }
            },
    A
            {
                @Override
                public AmberBase complement()
                {
                    return T;
                }
            },
    T
            {
                @Override
                public AmberBase complement()
                {
                    return A;
                }
            },
    C
            {
                @Override
                public AmberBase complement()
                {
                    return G;
                }
            },
    N
            {
                @Override
                public AmberBase complement()
                {
                    return N;
                }
            };

    public abstract AmberBase complement();
}
