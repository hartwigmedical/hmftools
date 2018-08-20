package com.hartwig.hmftools.strelka.mnv;

import java.util.Collection;
import java.util.List;
import java.util.Map;
import java.util.Optional;
import java.util.OptionalInt;
import java.util.function.Function;
import java.util.stream.Collectors;

import com.google.common.annotations.VisibleForTesting;
import com.google.common.collect.ArrayListMultimap;
import com.google.common.collect.Lists;
import com.google.common.collect.Multimap;

import org.apache.commons.lang3.tuple.Pair;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.immutables.value.Value;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.Genotype;
import htsjdk.variant.variantcontext.GenotypeBuilder;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.VariantContextBuilder;
import htsjdk.variant.vcf.VCFHeader;

@Value.Immutable
@Value.Style(allParameters = true,
             passAnnotations = { NotNull.class, Nullable.class })
public abstract class MNVMerger {
    private static final Logger LOGGER = LogManager.getLogger(MNVMerger.class);
    private static final String MNV_SET_VALUE = "mnvs";
    private static final String MNV_SET_KEY = "set";
    private static final String SOMATIC_PON_FIELD = "SOMATIC_PON_COUNT";
    private static final String GERMLINE_PON_FIELD = "GERMLINE_PON_COUNT";

    @NotNull
    abstract VCFHeader vcfHeader();

    @NotNull
    VariantContext mergeVariants(@NotNull final PotentialMNV mnv, @NotNull final Map<Integer, Character> gapReads) {
        final Map<Integer, Character> gapsForMnv =
                mnv.gapPositions().stream().collect(Collectors.toMap(Function.identity(), gapReads::get));
        return mergeVariants(mnv.variants(), gapsForMnv);
    }

    // MIVO: assumes variant list is sorted by start position and variants have only one sample (tumor)
    // MIVO: gaps will *ALWAYS* be added to the output variant (no checking is done here to make sure they are needed)
    @NotNull
    @VisibleForTesting
    VariantContext mergeVariants(@NotNull final List<VariantContext> variants, @NotNull final Map<Integer, Character> gapReads) {
        final List<Allele> alleles = createMnvAlleles(variants, gapReads);
        final VariantContext firstVariant = variants.get(0);
        final VariantContext lastVariant = variants.get(variants.size() - 1);
        final Map<String, Object> attributes = createMnvAttributes(variants);
        final String sampleName = firstVariant.getSampleNamesOrderedByName().get(0);
        final Genotype genotype = new GenotypeBuilder(sampleName, alleles).DP(mergeDP(variants)).AD(mergeAD(variants)).make();
        return new VariantContextBuilder(firstVariant.getSource(),
                firstVariant.getContig(),
                firstVariant.getStart(),
                lastVariant.getEnd(),
                alleles).genotypes(genotype).filters(firstVariant.getFilters()).attributes(attributes).make();
    }

    @NotNull
    private static List<Allele> createMnvAlleles(@NotNull final List<VariantContext> variants,
            @NotNull final Map<Integer, Character> gapReads) {
        final StringBuilder refBases = new StringBuilder();
        final StringBuilder altBases = new StringBuilder();
        int currentPosition = variants.get(0).getStart();
        for (final VariantContext variant : variants) {
            while (currentPosition != variant.getStart()) {
                final Character gapRead = gapReads.get(currentPosition);
                refBases.append(gapRead);
                altBases.append(gapRead);
                currentPosition++;
            }
            refBases.append(variant.getReference().getBaseString());
            altBases.append(variant.getAlternateAllele(0).getBaseString());
            currentPosition = variant.getStart() + variant.getReference().length();
        }
        final Allele ref = Allele.create(refBases.toString(), true);
        final Allele alt = Allele.create(altBases.toString(), false);
        return Lists.newArrayList(ref, alt);
    }

    private static int mergeDP(@NotNull final List<VariantContext> variants) {
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
    private static int[] mergeAD(@NotNull final List<VariantContext> variants) {
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
    private Map<String, Object> createMnvAttributes(@NotNull final List<VariantContext> variants) {
        final Map<String, Object> attributes;
        if (variants.stream().anyMatch(VariantContext::isIndel)) {
            attributes = mergeAttributes(variants.stream().filter(VariantContext::isIndel).collect(Collectors.toList()));
        } else {
            attributes = mergeAttributes(variants);
        }
        attributes.put(MNV_SET_KEY, MNV_SET_VALUE);
        return attributes;
    }

    @NotNull
    private Map<String, Object> mergeAttributes(@NotNull final List<VariantContext> variants) {
        final Multimap<String, Comparable> mergedAttributes = ArrayListMultimap.create();
        variants.stream()
                .flatMap(variant -> variant.fullyDecode(vcfHeader(), true).getAttributes().entrySet().stream())
                .filter(entry -> entry.getValue() instanceof Comparable)
                .forEach(entry -> mergedAttributes.put(entry.getKey(), (Comparable) entry.getValue()));
        return mergedAttributes.asMap().entrySet().stream().map(entry -> {
            final Optional<Comparable> minValue = determineFieldMin(entry, variants.size());
            return Pair.of(entry.getKey(), minValue);
        }).filter(pair -> pair.getValue().isPresent()).collect(Collectors.toMap(Pair::getKey, pair -> pair.getValue().get()));
    }

    private static Optional<Comparable> determineFieldMin(@NotNull final Map.Entry<String, Collection<Comparable>> entry,
            final int variantsSize) {
        if ((entry.getKey().equals(SOMATIC_PON_FIELD) || entry.getKey().equals(GERMLINE_PON_FIELD))
                && entry.getValue().size() != variantsSize) {
            return Optional.empty();
        }
        return entry.getValue().stream().sorted().findFirst();
    }
}
