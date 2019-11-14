package com.hartwig.hmftools.common.variant.strelka;

import java.util.List;
import java.util.function.Function;
import java.util.stream.Collectors;

import com.google.common.collect.Lists;

import org.jetbrains.annotations.NotNull;

import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.Genotype;
import htsjdk.variant.variantcontext.GenotypeBuilder;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.VariantContextBuilder;

public final class StrelkaAllelicDepth {

    private static final String TIR_FIELD = "TIR";
    private static final String TAR_FIELD = "TAR";
    private static final String SEPARATOR = ",";

    private StrelkaAllelicDepth() {
    }

    @NotNull
    public static VariantContext enrich(@NotNull final VariantContext context) {
        if (!context.isIndel() && !context.isSNP() ) {
            return context;
        }

        final VariantContextBuilder contextBuilder = new VariantContextBuilder(context).noGenotypes();

        final List<Allele> alleles = context.getAlleles();
        final Function<Allele, String> alleleKey = alleleKey(context);

        final List<Genotype> updatedGenotypes = Lists.newArrayList();
        for (Genotype genotype : context.getGenotypes()) {
            if (!genotype.hasAD() && hasRequiredAttributes(genotype, alleles, alleleKey)) {
                final GenotypeBuilder genotypeBuilder = new GenotypeBuilder(genotype).AD(readAD(genotype, alleles, alleleKey));
                updatedGenotypes.add(genotypeBuilder.make());
            } else {
                updatedGenotypes.add(genotype);
            }
        }

        return contextBuilder.genotypes(updatedGenotypes).make();
    }

    @NotNull
    public static String snpAlleleKey(@NotNull final Allele allele) {
        return allele.getBaseString() + "U";
    }

    @NotNull
    public static String indelAlleleKey(@NotNull final Allele allele) {
        return allele.isReference() ? TAR_FIELD : TIR_FIELD;
    }

    @NotNull
    public static int[] readAD(@NotNull final String sample, @NotNull final VariantContext variant) {
        final List<Allele> alleles = variant.getAlleles();
        final Function<Allele, String> alleleKey = alleleKey(variant);
        final Genotype tumorGenotype = variant.getGenotype(sample);
        return readAD(tumorGenotype, alleles, alleleKey);
    }

    @NotNull
    private static Function<Allele, String> alleleKey(@NotNull final VariantContext variant) {
        if (variant.isSNP()) {
            return StrelkaAllelicDepth::snpAlleleKey;
        } else if (variant.isIndel()) {
            return StrelkaAllelicDepth::indelAlleleKey;
        } else {
            throw new IllegalStateException("record is neither indel nor snp: " + variant);
        }
    }

    @NotNull
    private static int[] readAD(@NotNull final Genotype tumorGenotype, @NotNull final List<Allele> alleles,
            @NotNull Function<Allele, String> alleleKey) {
        final List<Integer> alleleAds = alleles.stream().map(allele -> {
            final String[] alleleADs = tumorGenotype.getExtendedAttribute(alleleKey.apply(allele), "0").toString().split(SEPARATOR);
            if (alleleADs.length > 0) {
                try {
                    return Integer.parseInt(alleleADs[0]);
                } catch (final NumberFormatException e) {
                    return 0;
                }
            }
            return 0;
        }).collect(Collectors.toList());
        return alleleAds.stream().mapToInt(Integer::intValue).toArray();
    }

    private static boolean hasRequiredAttributes(@NotNull final Genotype tumorGenotype, @NotNull final List<Allele> alleles,
            @NotNull Function<Allele, String> alleleKey) {
        return alleles.stream().allMatch(x -> tumorGenotype.hasExtendedAttribute(alleleKey.apply(x)));
    }
}
