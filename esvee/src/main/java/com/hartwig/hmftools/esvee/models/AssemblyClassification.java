package com.hartwig.hmftools.esvee.models;

import static com.hartwig.hmftools.common.sv.StructuralVariantType.BND;
import static com.hartwig.hmftools.common.sv.StructuralVariantType.INV;
import static com.hartwig.hmftools.common.sv.StructuralVariantType.SGL;

import com.hartwig.hmftools.common.sv.StructuralVariantType;

public class AssemblyClassification
{
    public final StructuralVariantType Type;

    /** The length of this event, if applicable. 0 otherwise */
    public final int Length;

    public AssemblyClassification(final StructuralVariantType type, final int length)
    {
        Type = type;
        Length = length;
    }

    @Override
    public String toString()
    {
        if(Type != SGL && Type != BND && Type != INV)
            return Type.toString();
        else
            return Type.name() + " " + Length;
    }
}
