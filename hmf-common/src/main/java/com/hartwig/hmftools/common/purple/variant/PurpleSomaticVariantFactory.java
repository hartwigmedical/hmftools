package com.hartwig.hmftools.common.purple.variant;

import static htsjdk.tribble.AbstractFeatureReader.getFeatureReader;

import java.io.IOException;
import java.util.List;
import java.util.stream.Collectors;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.variant.VariantType;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.jetbrains.annotations.NotNull;

import htsjdk.tribble.AbstractFeatureReader;
import htsjdk.tribble.readers.LineIterator;
import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.Genotype;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.filter.CompoundFilter;
import htsjdk.variant.variantcontext.filter.PassingVariantFilter;
import htsjdk.variant.variantcontext.filter.SnpFilter;
import htsjdk.variant.vcf.VCFCodec;
import htsjdk.variant.vcf.VCFHeader;

public class PurpleSomaticVariantFactory {

    private static final Logger LOGGER = LogManager.getLogger(PurpleSomaticVariantFactory.class);

    private final CompoundFilter filter;

    public PurpleSomaticVariantFactory() {
        final CompoundFilter filter = new CompoundFilter(true);
        filter.add(new PassingVariantFilter());
        filter.add(new SnpFilter());

        this.filter = filter;
    }

    public List<PurpleSomaticVariant> fromVCFFile(@NotNull final String sample, @NotNull final String vcfFile) throws IOException {
        final List<PurpleSomaticVariant> variants = Lists.newArrayList();

        try (final AbstractFeatureReader<VariantContext, LineIterator> reader = getFeatureReader(vcfFile, new VCFCodec(), false)) {

            if (!sampleInFile(sample, (VCFHeader) reader.getHeader())) {
                throw new IllegalArgumentException("Sample " + sample + " not found in vcf file " + vcfFile);
            }

            for (final VariantContext context : reader.iterator()) {
                if (filter.test(context)) {
                    final Genotype genotype = context.getGenotype(sample);
                    if (genotype.hasAD()) {
                        variants.add(ImmutablePurpleSomaticVariant.builder()
                                .chromosome(context.getContig())
                                .position(context.getStart())
                                .ref(context.getReference().getBaseString())
                                .alt(alt(context))
                                .alleleFrequency(allelicFrequency(genotype.getAD()))
                                .filter("PASS")
                                .type(VariantType.SNP)
                                .build());
                    } else {
                        LOGGER.warn("Missing AD flag.");
                    }

                }
            }
        }

        return variants;
    }

    private boolean sampleInFile(@NotNull final String sample, @NotNull final VCFHeader header) {
        return header.getSampleNamesInOrder().stream().anyMatch(x -> x.equals(sample));
    }

    @NotNull
    private static String alt(@NotNull final VariantContext context) {
        return String.join(",", context.getAlternateAlleles().stream().map(Allele::toString).collect(Collectors.toList()));
    }

    private static double allelicFrequency(int[] afFields) {

        int totalReadCount = 0;
        final int alleleReadCount = afFields[1];
        for (final int afField : afFields) {
            totalReadCount += afField;
        }
        return alleleReadCount / (double) totalReadCount;
    }
}
