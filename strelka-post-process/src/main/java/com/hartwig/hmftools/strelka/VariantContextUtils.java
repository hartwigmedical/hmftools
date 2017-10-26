package com.hartwig.hmftools.strelka;

import java.util.List;

import com.google.common.collect.Lists;

import org.jetbrains.annotations.NotNull;

import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.Genotype;
import htsjdk.variant.variantcontext.GenotypeBuilder;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.VariantContextBuilder;

public final class VariantContextUtils {
    private VariantContextUtils() {
    }

    @NotNull
    public static List<VariantContext> splitMultiAlleleVariant(@NotNull final VariantContext variant) {
        final List<VariantContext> variants = Lists.newArrayList();
        final int[] ad = variant.getGenotype(0).getAD();
        final int dp = variant.getGenotype(0).getDP();
        for (int index = 0; index < variant.getAlternateAlleles().size(); index++) {
            final Allele alt = variant.getAlternateAllele(index);
            final List<Allele> alleles = Lists.newArrayList(variant.getReference(), alt);
            final Genotype genotype = new GenotypeBuilder(variant.getSampleNamesOrderedByName().get(0), alleles).DP(dp)
                    .AD(new int[] { ad[0], ad[index + 1] })
                    .make();
            variants.add(new VariantContextBuilder(variant).alleles(alleles).genotypes(genotype).make());
        }
        return variants;
    }
}
