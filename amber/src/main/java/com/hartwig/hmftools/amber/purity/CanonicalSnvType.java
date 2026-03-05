package com.hartwig.hmftools.amber.purity;

import static com.hartwig.hmftools.common.amber.AmberBase.A;
import static com.hartwig.hmftools.common.amber.AmberBase.C;
import static com.hartwig.hmftools.common.amber.AmberBase.G;
import static com.hartwig.hmftools.common.amber.AmberBase.T;

import com.google.common.base.Preconditions;
import com.hartwig.hmftools.common.amber.AmberBase;

public enum CanonicalSnvType
{
    C_A,
    C_G,
    C_T,
    T_A,
    T_C,
    T_G;

    public static CanonicalSnvType type(AmberBase ref, AmberBase alt)
    {
        Preconditions.checkArgument(ref != alt);
        Preconditions.checkArgument(ref != AmberBase.N);
        Preconditions.checkArgument(alt != AmberBase.N);
        AmberBase canonicalRef = ref;
        AmberBase canonicalAlt = alt;
        if(ref == A || ref == G)
        {
            canonicalRef = ref.complement();
            canonicalAlt = alt.complement();
        }

        if(canonicalRef == C)
        {
            switch(canonicalAlt)
            {
                case A:
                    return C_A;
                case G:
                    return C_G;
                case T:
                    return C_T;
            }
        }
        if(canonicalRef == T)
        {
            switch(canonicalAlt)
            {
                case A:
                    return T_A;
                case C:
                    return T_C;
                case G:
                    return T_G;
            }
        }
        return null;
    }
}
