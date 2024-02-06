package com.hartwig.hmftools.wisp.purity.variant;

import com.hartwig.hmftools.common.variant.VariantType;

public class ProbeVariant
{
    public final String Chromosome;
    public final int Position;
    public final String Ref;
    public final String Alt;
    public final VariantType Type;

    public ProbeVariant(final String chromosome, final int position, final String ref, final String alt, final VariantType type)
    {
        Chromosome = chromosome;
        Position = position;
        Ref = ref;
        Alt = alt;
        Type = type;
    }

    public boolean matches(final SomaticVariant somaticVariant)
    {
        return somaticVariant.Chromosome.equals(Chromosome) && somaticVariant.Position == Position
                && somaticVariant.Ref.equals(Ref) && somaticVariant.Alt.equals(Alt) && somaticVariant.Type == Type;
    }
}
