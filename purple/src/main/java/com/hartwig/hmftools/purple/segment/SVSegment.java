package com.hartwig.hmftools.purple.segment;

import com.hartwig.hmftools.common.genome.position.GenomePosition;
import com.hartwig.hmftools.common.sv.StructuralVariantType;

public class SVSegment implements GenomePosition
{
    public final String Chromosome;
    public final int Position;
    public final StructuralVariantType Type;

    public SVSegment(final String chromosome, final int position, final StructuralVariantType type)
    {
        Chromosome = chromosome;
        Position = position;
        Type = type;
    }

    @Override
    public String chromosome() { return Chromosome; }

    @Override
    public int position()
    {
        return Position;
    }
}
