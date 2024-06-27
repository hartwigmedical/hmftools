package com.hartwig.hmftools.common.variant;

import static java.lang.String.format;

import com.google.common.base.Preconditions;

import org.jetbrains.annotations.NotNull;

import htsjdk.variant.variantcontext.Genotype;

public class AllelicDepth
{
    public final int TotalReadCount;
    public final int AlleleReadCount;

    public AllelicDepth(final int totalReadCount, final int alleleReadCount)
    {
        TotalReadCount = totalReadCount;
        AlleleReadCount = alleleReadCount;
    }

    public static final AllelicDepth NO_DEPTH = new AllelicDepth(0, 0);

    public double alleleFrequency()
    {
        return TotalReadCount > 0 ? AlleleReadCount / (double)TotalReadCount : 0;
    }

    public String toString() { return format("%d/%d", AlleleReadCount, TotalReadCount); }

    public static boolean containsAllelicDepth(final Genotype genotype)
    {
        return genotype != null && genotype.hasAD() && genotype.getAD().length > 1;
    }

    public static AllelicDepth fromGenotype(final Genotype genotype)
    {
        Preconditions.checkArgument(genotype.hasAD());
        int[] adFields = genotype.getAD();
        final int alleleReadCount = adFields[1];
        int totalReadCount = totalReadCount(genotype);
        return new AllelicDepth(totalReadCount, alleleReadCount);
    }

    public static int totalReadCount(final Genotype genotype)
    {
        // Note: this is a workaround of strelka's DP being only Tier 1
        return genotype.hasDP() ? Math.max(genotype.getDP(), sumReadCount(genotype.getAD())) : sumReadCount(genotype.getAD());
    }

    public static int sumReadCount(int[] adFields)
    {
        int totalReadCount = 0;
        for(int afField : adFields)
        {
            totalReadCount += afField;
        }
        return totalReadCount;
    }
}
