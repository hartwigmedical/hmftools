package com.hartwig.hmftools.strelka;

import java.util.List;
import java.util.Map;
import java.util.OptionalInt;
import java.util.stream.Collectors;

import com.google.common.collect.Lists;
import com.google.common.collect.Maps;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.jetbrains.annotations.NotNull;

import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.Genotype;
import htsjdk.variant.variantcontext.GenotypeBuilder;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.VariantContextBuilder;

public class MNVMerger {
    private static final Logger LOGGER = LogManager.getLogger(MNVMerger.class);
    private static final String SET_VALUE = "mnvs";
    private static final String SET_KEY = "set";

    private MNVMerger() {
    }

    @NotNull
    public static VariantContext mergeVariants(@NotNull final String sampleName, @NotNull final List<VariantContext> variants,
            @NotNull final Map<Integer, Character> gapReads) {
        //MIVO: assumes variant list is sorted by start position
        final List<Allele> alleles = createMnvAlleles(variants, gapReads);
        final Genotype genotype = new GenotypeBuilder(sampleName, alleles).DP(mergeDP(variants)).AD(mergeAD(variants)).make();
        final VariantContext firstVariant = variants.get(0);
        final VariantContext lastVariant = variants.get(variants.size() - 1);
        final Map<String, Object> attributes = createMnvAttributes(variants);
        return new VariantContextBuilder("", firstVariant.getContig(), firstVariant.getStart(), lastVariant.getEnd(), alleles).genotypes(
                genotype).filters(firstVariant.getFilters()).attributes(attributes).make();
    }

    @NotNull
    static List<Allele> createMnvAlleles(@NotNull final List<VariantContext> variants, @NotNull final Map<Integer, Character> gapReads) {
        final StringBuilder refBases = new StringBuilder();
        final StringBuilder altBases = new StringBuilder();
        variants.forEach(variant -> {
            final int positionBeforeVariant = variant.getStart() - 1;
            if (gapReads.containsKey(positionBeforeVariant)) {
                final Character gapRead = gapReads.get(positionBeforeVariant);
                refBases.append(gapRead);
                altBases.append(gapRead);
            }
            refBases.append(variant.getReference().getBaseString());
            altBases.append(variant.getAlternateAllele(0).getBaseString());
        });
        final Allele ref = Allele.create(refBases.toString(), true);
        final Allele alt = Allele.create(altBases.toString(), false);
        return Lists.newArrayList(ref, alt);
    }

    @NotNull
    static int mergeDP(@NotNull final List<VariantContext> variants) {
        final OptionalInt min = variants.stream().mapToInt(variant -> variant.getGenotype(0).getDP()).min();
        if (min.isPresent()) {
            return min.getAsInt();
        } else {
            LOGGER.warn("No min DP found for variants: {}",
                    variants.stream().map(variant -> variant.getContig() + ":" + variant.getStart()).collect(Collectors.joining(",")));
            return 0;
        }
    }

    @NotNull
    static int[] mergeAD(@NotNull final List<VariantContext> variants) {
        int[] ads = new int[] { Integer.MAX_VALUE, Integer.MAX_VALUE };
        variants.forEach(variant -> {
            int[] variantADs = variant.getGenotype(0).getAD();
            if (variantADs[0] < ads[0]) {
                ads[0] = variantADs[0];
            }
            if (variantADs[1] < ads[1]) {
                ads[1] = variantADs[1];
            }
        });
        return ads;
    }

    @NotNull
    static Map<String, Object> createMnvAttributes(@NotNull final List<VariantContext> variants) {
        final Map<String, Object> attributes;
        if (variants.stream().anyMatch(VariantContext::isIndel)) {
            attributes = mergeAttributes(variants.stream().filter(VariantContext::isIndel).collect(Collectors.toList()));
        } else {
            attributes = mergeAttributes(variants);
        }
        attributes.put(SET_KEY, SET_VALUE);
        return attributes;
    }

    @NotNull
    static Map<String, Object> mergeAttributes(@NotNull final List<VariantContext> variants) {
        final Map<String, Object> mergedAttributes = Maps.newHashMap();
        variants.forEach(variant -> variant.getAttributes().forEach((key, value) -> {
            if (mergedAttributes.containsKey(key)) {
                if (value instanceof Comparable && mergedAttributes.get(key) instanceof Comparable) {
                    final Comparable entryValue = (Comparable) value;
                    final Comparable mergedValue = (Comparable) mergedAttributes.get(key);
                    if (entryValue.compareTo(mergedValue) < 1) {
                        mergedAttributes.put(key, value);
                    }
                }
            } else {
                mergedAttributes.put(key, value);
            }
        }));
        return mergedAttributes;
    }
}
