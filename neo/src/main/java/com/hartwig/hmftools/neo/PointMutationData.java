package com.hartwig.hmftools.neo;

import com.hartwig.hmftools.common.variant.CodingEffect;

public class PointMutationData
{
    public final String Chromosome;
    public final int Position;
    public final String Ref;
    public final String Alt;
    public final String Gene;
    public final CodingEffect Effect;
    public final double CopyNumber;
    public final int LocalPhaseSet;

    public PointMutationData(
            final String chromosome, final int position, final String ref, final String alt, final String gene,
            final CodingEffect effect, double copyNumber, int localPhaseSet)
    {
        Chromosome = chromosome;
        Position = position;
        Ref = ref;
        Alt = alt;
        Gene = gene;
        Effect = effect;
        CopyNumber = copyNumber;
        LocalPhaseSet = localPhaseSet;
    }
}
