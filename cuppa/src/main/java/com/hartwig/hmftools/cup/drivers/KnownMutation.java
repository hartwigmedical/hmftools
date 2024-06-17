package com.hartwig.hmftools.cup.drivers;

import static com.hartwig.hmftools.common.region.BaseRegion.positionWithin;
import static com.hartwig.hmftools.common.utils.sv.StartEndIterator.SE_END;
import static com.hartwig.hmftools.common.utils.sv.StartEndIterator.SE_START;

import com.hartwig.hmftools.common.variant.VariantType;

public class KnownMutation
{
    public final String Gene;
    public final VariantType Type;
    public final String Ref;
    public final String Alt;
    public final int[] Position;

    public KnownMutation(final String gene, final VariantType type, final String ref, final String alt, final int posStart, final int posEnd)
    {
        Gene = gene;
        Type = type;
        Ref = ref;
        Alt = alt;
        Position = new int[] { posStart, posEnd };
    }

    public boolean matches(final String gene, final VariantType type, final String ref, final String alt, final int position)
    {
        return Gene.equals(gene) && Type == type
            && (Ref.isEmpty() || Ref.equals(ref)) && (Alt.isEmpty() || Alt.equals(alt))
            && positionWithin(position, Position[SE_START], Position[SE_END]);
    }
}
