package com.hartwig.hmftools.common.variant;

import static htsjdk.tribble.AbstractFeatureReader.getFeatureReader;
import static htsjdk.variant.variantcontext.Genotype.PHASED_ALLELE_SEPARATOR;
import static htsjdk.variant.variantcontext.Genotype.UNPHASED_ALLELE_SEPARATOR;

import java.io.IOException;
import java.util.List;
import java.util.StringJoiner;
import java.util.function.Predicate;
import java.util.stream.Collectors;

import com.google.common.collect.Lists;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.jetbrains.annotations.NotNull;

import htsjdk.tribble.AbstractFeatureReader;
import htsjdk.tribble.readers.LineIterator;
import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.Genotype;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.vcf.VCFCodec;
import htsjdk.variant.vcf.VCFHeader;

public class GermlineVariantFactory2 {
    private static final Logger LOGGER = LogManager.getLogger(GermlineVariantFactory2.class);

    private Predicate<VariantContext> variantContextFilter = x -> true;
    private Predicate<Genotype> normalGenotypeFilter = this::standardGenotypeFilter;
    private Predicate<Genotype> tumorGenotypeFilter = this::standardGenotypeFilter;

    public void addVariantFilter(Predicate<VariantContext> filter) {
        variantContextFilter = variantContextFilter.and(filter);
    }

    public void addNormalGenotypeFilter(Predicate<Genotype> filter) {
        normalGenotypeFilter = normalGenotypeFilter.and(filter);
    }

    public void addTumorGenotypeFilter(Predicate<Genotype> filter) {
        tumorGenotypeFilter = tumorGenotypeFilter.and(filter);
    }

    public List<GermlineVariant> fromVCFFile(@NotNull String vcfFile) throws IOException {
        final List<GermlineVariant> variants = Lists.newArrayList();

        try (final AbstractFeatureReader<VariantContext, LineIterator> reader = getFeatureReader(vcfFile, new VCFCodec(), false)) {

            VCFHeader header = (VCFHeader) reader.getHeader();

            final String normalSample = header.getSampleNamesInOrder().get(0);
            final String tumorSample = header.getSampleNamesInOrder().get(1);

            for (final VariantContext context : reader.iterator()) {

                if (variantContextFilter.test(context)) {
                    final Genotype normalGenotype = context.getGenotype(normalSample);
                    final Genotype tumorGenotype = context.getGenotype(tumorSample);
                    if (normalGenotypeFilter.test(normalGenotype) && tumorGenotypeFilter.test(tumorGenotype)) {
                        final String chromosome = context.getContig();
                        final long start = context.getStart();

                        variants.add(ImmutableGermlineVariant.builder()
                                .chromosome(chromosome)
                                .position(start)
                                .ref(context.getReference().getBaseString())
                                .alt(String.join(",",
                                        context.getAlternateAlleles().stream().map(Allele::toString).collect(Collectors.toList())))
                                .refData(create(context, normalGenotype))
                                .tumorData(create(context, tumorGenotype))
                                .type(variantType(context.getType()))
                                .filter("PASS") //TODO: Fix this
                                .build());
                    }
                }
            }
        }

        return variants;
    }

    private static VariantType variantType(VariantContext.Type type) {
        switch (type) {
            case SNP:
                return VariantType.SNP;
            case INDEL:
                return VariantType.INDEL;
            default:
                return VariantType.UNDEFINED;
        }
    }

    private static GermlineSampleData create(VariantContext context, Genotype genotype) {
        final AlleleFrequencyData frequencyData = createAlleleFrequencyData(genotype.getAD());
        return ImmutableGermlineSampleData.builder()
                .genoType(getGenotypeString(context, genotype))
                .totalReadCount(frequencyData.totalReadCount())
                .alleleReadCount(frequencyData.alleleReadCount())
                .combinedDepth(genotype.getDP())
                .build();

    }

    private static String getGenotypeString(VariantContext context, Genotype genotype) {
        final String separator = genotype.isPhased() ? PHASED_ALLELE_SEPARATOR : UNPHASED_ALLELE_SEPARATOR;
        StringJoiner joiner = new StringJoiner(separator);
        for (Allele genoTypeAllele : genotype.getAlleles()) {
            int index = context.getAlleleIndex(genoTypeAllele);
            joiner.add("" + (index == -1 ? "." : index));
        }

        return joiner.toString();
    }

    private static AlleleFrequencyData createAlleleFrequencyData(int[] afFields) {

        int totalReadCount = 0;
        final int alleleReadCount = afFields[1];
        for (final int afField : afFields) {
            totalReadCount += afField;
        }
        return new AlleleFrequencyData(alleleReadCount, totalReadCount);
    }

    private boolean standardGenotypeFilter(Genotype genotype) {
        if (!genotype.hasAD()) {
            LOGGER.warn("Could not parse allele frequencies for genotype: " + genotype);
            return false;
        }

        if (!genotype.hasDP()) {
            LOGGER.warn("Could not parse combined depth for genotype: " + genotype);
            return false;
        }

        return true;
    }
}
