package com.hartwig.hmftools.common.variant;

import static org.junit.Assert.assertEquals;

import org.junit.Test;

import htsjdk.variant.variantcontext.Genotype;
import htsjdk.variant.variantcontext.GenotypeBuilder;

public class AllelicDepthTest
{
    @Test
    public void allelicDepthFromAd()
    {
        final Genotype genotype = new GenotypeBuilder("SAMPLE").AD(new int[] { 6, 4 }).make();
        final AllelicDepth victim = AllelicDepth.fromGenotype(genotype);
        assertEquals(10, victim.TotalReadCount);
        assertEquals(0.4, victim.alleleFrequency(), 0.01);
    }

    @Test
    public void allelicDepthFromDpIfAvailable()
    {
        final Genotype genotype = new GenotypeBuilder("SAMPLE").AD(new int[] { 6, 4 }).DP(20).make();
        final AllelicDepth victim = AllelicDepth.fromGenotype(genotype);
        assertEquals(20, victim.TotalReadCount);
        assertEquals(0.2, victim.alleleFrequency(), 0.01);
    }

    @Test
    public void useADIfLargerThanDP()
    {
        final Genotype genotype = new GenotypeBuilder("SAMPLE").AD(new int[] { 6, 4 }).DP(5).make();
        final AllelicDepth victim = AllelicDepth.fromGenotype(genotype);
        assertEquals(10, victim.TotalReadCount);
        assertEquals(0.4, victim.alleleFrequency(), 0.01);
    }
}
