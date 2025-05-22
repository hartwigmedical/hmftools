package com.hartwig.hmftools.esvee.vcfcompare;

import static com.hartwig.hmftools.common.region.BaseRegion.positionsOverlap;
import static com.hartwig.hmftools.common.sv.StructuralVariantType.DUP;
import static com.hartwig.hmftools.common.sv.StructuralVariantType.INS;

import com.hartwig.hmftools.common.sv.StructuralVariantType;

public enum CoordMatchType
{
    EXACT,
    DUP_INS,
    HOMOLOGY,
    ONE_SIDE,
    SV_AND_SGL,
    NONE;

    public static CoordMatchType determineMatchType(final Breakend first, final Breakend second)
    {
        CoordMatchType startMatch = NONE;

        if(exactMatch(first, second))
        {
            startMatch = EXACT;
        }
        else if(homologyMatch(first, second))
        {
            startMatch = HOMOLOGY;
        }

        if(first.isSgl() && second.isSgl())
            return startMatch;

        if(first.isSgl() != second.isSgl())
            return startMatch == NONE ? NONE : SV_AND_SGL;

        CoordMatchType endMatch = NONE;

        if(exactMatch(first.otherBreakend(), second.otherBreakend()))
        {
            endMatch = EXACT;
        }
        else if(homologyMatch(first.otherBreakend(), second.otherBreakend()))
        {
            endMatch = HOMOLOGY;
        }

        if(startMatch == NONE && endMatch == NONE)
            return NONE;

        if(isDupInsPair(first.type(), second.type()))
        {
            if(startMatch == EXACT || endMatch == EXACT)
                return DUP_INS;
        }

        if(startMatch == NONE || endMatch == NONE)
            return ONE_SIDE;

        if(startMatch == EXACT && endMatch == EXACT)
            return EXACT;

        return HOMOLOGY;
    }

    private static boolean isDupInsPair(final StructuralVariantType first, final StructuralVariantType second)
    {
        return (first == DUP && second == INS) || (first == INS && second == DUP);
    }

    public static boolean exactMatch(final Breakend first, final Breakend second)
    {
        return first.Chromosome.equals(second.Chromosome) && first.Position == second.Position && first.Orient == second.Orient;
    }

    public static boolean homologyMatch(final Breakend first, final Breakend second)
    {
        if(!first.Chromosome.equals(second.Chromosome) || first.Orient != second.Orient)
            return false;

        return positionsOverlap(first.minPosition(), first.maxPosition(), second.minPosition(), second.maxPosition());
    }
}
