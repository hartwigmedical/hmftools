package com.hartwig.hmftools.svtools.neo;

public class PointMutationData
{
    public final String Chromosome;
    public final int Position;
    public final String Ref;
    public final String Alt;
    public final String Gene;
    public final int LocalPhaseSet;

    public PointMutationData(
            final String chromosome, final int position, final String ref, final String alt, final String gene, int localPhaseSet)
    {
        Chromosome = chromosome;
        Position = position;
        Ref = ref;
        Alt = alt;
        Gene = gene;
        LocalPhaseSet = localPhaseSet;
    }
}
