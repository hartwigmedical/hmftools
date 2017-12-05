package com.hartwig.hmftools.common.purple.variant;

import static htsjdk.tribble.AbstractFeatureReader.getFeatureReader;

import java.io.IOException;
import java.util.List;
import java.util.stream.Collectors;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.variant.AlleleFrequencyData;
import com.hartwig.hmftools.common.variant.ChromosomeFilter;
import com.hartwig.hmftools.common.variant.VariantFactoryFunctions;
import com.hartwig.hmftools.common.variant.VariantType;

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

    @NotNull
    private final CompoundFilter filter;

    public PurpleSomaticVariantFactory() {
        final CompoundFilter filter = new CompoundFilter(true);
        filter.add(new PassingVariantFilter());
        filter.add(new SnpFilter());
        filter.add(new ChromosomeFilter());
        filter.add(PurpleSomaticVariantFactory::ntFilter);
        filter.add(PurpleSomaticVariantFactory::sgtFilter);

        this.filter = filter;
    }

    @NotNull
    public List<PurpleSomaticVariant> fromVCFFile(@NotNull final String sample, @NotNull final String vcfFile) throws IOException {
        final List<PurpleSomaticVariant> variants = Lists.newArrayList();

        try (final AbstractFeatureReader<VariantContext, LineIterator> reader = getFeatureReader(vcfFile, new VCFCodec(), false)) {

            final VCFHeader header = (VCFHeader) reader.getHeader();
            if (!sampleInFile(sample, header)) {
                throw new IllegalArgumentException("Sample " + sample + " not found in vcf file " + vcfFile);
            }

            if (!header.hasFormatLine("AD")) {
                throw new IllegalArgumentException("Allelic depths is a required format field in vcf file " + vcfFile);
            }

            for (final VariantContext context : reader.iterator()) {
                if (filter.test(context)) {
                    final Genotype genotype = context.getGenotype(sample);
                    if (genotype.hasAD()) {
                        final AlleleFrequencyData frequencyData = VariantFactoryFunctions.determineAlleleFrequencies(genotype);
                        variants.add(ImmutablePurpleSomaticVariant.builder()
                                .chromosome(context.getContig())
                                .position(context.getStart())
                                .ref(context.getReference().getBaseString())
                                .alt(alt(context))
                                .alleleReadCount(frequencyData.alleleReadCount())
                                .totalReadCount(frequencyData.totalReadCount())
                                .filter("PASS")
                                .type(VariantType.SNP)
                                .build());
                    }
                }
            }
        }

        return variants;
    }

    private static boolean sampleInFile(@NotNull final String sample, @NotNull final VCFHeader header) {
        return header.getSampleNamesInOrder().stream().anyMatch(x -> x.equals(sample));
    }

    @NotNull
    private static String alt(@NotNull final VariantContext context) {
        return String.join(",", context.getAlternateAlleles().stream().map(Allele::toString).collect(Collectors.toList()));
    }

    private static boolean ntFilter(@NotNull final VariantContext record) {
        return !record.getCommonInfo().hasAttribute("NT") || record.getCommonInfo().getAttribute("NT").equals("ref");
    }

    private static boolean sgtFilter(@NotNull final VariantContext record) {
        if (record.getCommonInfo().hasAttribute("SGT")) {
            final String[] sgt = record.getCommonInfo().getAttributeAsString("SGT", "GC->AT").split("->");
            if (sgt.length == 2 && sgt[0].equals(sgt[1])) {
                return false;
            }
        }
        return true;
    }
}
