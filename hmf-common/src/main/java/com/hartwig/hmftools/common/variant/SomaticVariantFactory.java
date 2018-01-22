package com.hartwig.hmftools.common.variant;

import static htsjdk.tribble.AbstractFeatureReader.getFeatureReader;

import java.io.IOException;
import java.util.Arrays;
import java.util.Collections;
import java.util.List;
import java.util.Optional;
import java.util.StringJoiner;
import java.util.stream.Collectors;

import com.google.common.annotations.VisibleForTesting;
import com.google.common.base.Preconditions;
import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.io.path.PathExtensionFinder;
import com.hartwig.hmftools.common.variant.filter.ChromosomeFilter;

import org.jetbrains.annotations.NotNull;

import htsjdk.tribble.AbstractFeatureReader;
import htsjdk.tribble.readers.LineIterator;
import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.Genotype;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.filter.CompoundFilter;
import htsjdk.variant.variantcontext.filter.VariantContextFilter;
import htsjdk.variant.vcf.VCFCodec;
import htsjdk.variant.vcf.VCFHeader;

public class SomaticVariantFactory {
    private static final String DBSNP_IDENTIFIER = "rs";
    private static final String COSMIC_IDENTIFIER = "COSM";
    private static final String ID_SEPARATOR = ";";
    private static final String MAPPABILITY_TAG = "MAPPABILITY";

    private static final String CALLER_ALGO_IDENTIFIER = "set";
    private static final String CALLER_ALGO_SEPARATOR = "-";
    private static final String CALLER_FILTERED_IDENTIFIER = "filterIn";
    private static final String CALLER_INTERSECTION_IDENTIFIER = "Intersection";

    @NotNull
    private final CompoundFilter filter;

    public SomaticVariantFactory() {
        this(new VariantContextFilter[0]);
    }

    public SomaticVariantFactory(VariantContextFilter... filters) {
        final CompoundFilter filter = new CompoundFilter(true);
        filter.add(new ChromosomeFilter());
        filter.addAll(Arrays.asList(filters));
        this.filter = filter;
    }

    @NotNull
    public List<SomaticVariant> fromVCFFile(@NotNull final String sample, @NotNull final String basePath,
            @NotNull final String fileExtension) throws IOException {
        return fromVCFFile(sample, PathExtensionFinder.build().findPath(basePath, fileExtension).toString());
    }

    @NotNull
    public List<SomaticVariant> fromVCFFile(@NotNull final String sample, @NotNull final String vcfFile) throws IOException {
        final List<SomaticVariant> variants = Lists.newArrayList();

        try (final AbstractFeatureReader<VariantContext, LineIterator> reader = getFeatureReader(vcfFile, new VCFCodec(), false)) {

            final VCFHeader header = (VCFHeader) reader.getHeader();
            if (!sampleInFile(sample, header)) {
                throw new IllegalArgumentException("Sample " + sample + " not found in vcf file " + vcfFile);
            }

            if (!header.hasFormatLine("AD")) {
                throw new IllegalArgumentException("Allelic depths is a required format field in vcf file " + vcfFile);
            }

            for (final VariantContext context : reader.iterator()) {
                createVariant(sample, context).ifPresent(variants::add);
            }
        }

        return variants;
    }

    @VisibleForTesting
    @NotNull
    Optional<SomaticVariant> createVariant(@NotNull final String sample, @NotNull final VariantContext context) {
        if (filter.test(context)) {
            final Genotype genotype = context.getGenotype(sample);
            if (genotype.hasAD() && genotype.getAD().length > 1) {
                final AlleleFrequencyData frequencyData = determineAlleleFrequencies(genotype);
                if (frequencyData.totalReadCount() > 0) {
                    SomaticVariantImpl.Builder builder = new SomaticVariantImpl.Builder().chromosome(context.getContig())
                            .annotations(Collections.emptyList())
                            .position(context.getStart())
                            .ref(context.getReference().getBaseString())
                            .alt(alt(context))
                            .alleleReadCount(frequencyData.alleleReadCount())
                            .totalReadCount(frequencyData.totalReadCount())
                            .totalReadCount(frequencyData.totalReadCount())
                            .mappability(context.getAttributeAsDouble(MAPPABILITY_TAG, 0));

                    attachCallers(builder, context);
                    attachAnnotations(builder, context);
                    attachFilter(builder, context);
                    attachID(builder, context);
                    attachType(builder, context);
                    return Optional.of(builder.build());
                }
            }
        }
        return Optional.empty();
    }

    private static void attachAnnotations(@NotNull final SomaticVariantImpl.Builder builder,
            @NotNull VariantContext context) {
        builder.annotations(VariantAnnotationFactory.fromContext(context));
    }

    private static void attachCallers(@NotNull final SomaticVariantImpl.Builder builder,
            @NotNull VariantContext context) {
        if (context.getCommonInfo().hasAttribute(CALLER_ALGO_IDENTIFIER)) {
            builder.callers(extractCallers(context.getCommonInfo().getAttributeAsString(CALLER_ALGO_IDENTIFIER, "")));
        } else {
            builder.callers(Collections.emptyList());
        }
    }

    @NotNull
    private static List<String> extractCallers(@NotNull final String callers) {
        final String[] allCallers = callers.split(CALLER_ALGO_SEPARATOR);
        final List<String> finalCallers = Lists.newArrayList();
        if (allCallers.length > 0 && allCallers[0].equals(CALLER_INTERSECTION_IDENTIFIER)) {
            finalCallers.addAll(SomaticVariantConstants.ALL_CALLERS);
        } else {
            finalCallers.addAll(Arrays.stream(allCallers)
                    .filter(caller -> !caller.startsWith(CALLER_FILTERED_IDENTIFIER))
                    .collect(Collectors.toList()));
        }
        return finalCallers;
    }

    private static void attachFilter(@NotNull final SomaticVariantImpl.Builder builder,
            @NotNull VariantContext context) {
        if (context.isFiltered()) {
            StringJoiner joiner = new StringJoiner(";");
            context.getFilters().forEach(joiner::add);
            builder.filter(joiner.toString());
        }

        builder.filter("PASS");
    }

    private static void attachType(@NotNull final SomaticVariantImpl.Builder builder,
            @NotNull VariantContext context) {
        switch (context.getType()) {
            case MNP:
                builder.type(VariantType.MNP);
                break;
            case SNP:
                builder.type(VariantType.SNP);
                break;
            case INDEL:
                builder.type(VariantType.INDEL);
                break;
            default:
                builder.type(VariantType.UNDEFINED);
        }
    }

    private static void attachID(@NotNull SomaticVariantImpl.Builder builder, @NotNull VariantContext context) {
        final String ID = context.getID();
        if (!ID.isEmpty()) {
            final String[] ids = ID.split(ID_SEPARATOR);
            for (final String id : ids) {
                if (id.contains(DBSNP_IDENTIFIER)) {
                    builder.dbsnpID(id);
                } else if (id.contains(COSMIC_IDENTIFIER)) {
                    builder.cosmicID(id);
                }
            }
        }
    }

    private static boolean sampleInFile(@NotNull final String sample, @NotNull final VCFHeader header) {
        return header.getSampleNamesInOrder().stream().anyMatch(x -> x.equals(sample));
    }

    @NotNull
    private static String alt(@NotNull final VariantContext context) {
        return String.join(",", context.getAlternateAlleles().stream().map(Allele::toString).collect(Collectors.toList()));
    }

    @NotNull
    private static AlleleFrequencyData determineAlleleFrequencies(@NotNull final Genotype genotype) {
        Preconditions.checkArgument(genotype.hasAD());

        int[] adFields = genotype.getAD();
        int totalReadCount = 0;
        final int alleleReadCount = adFields[1];
        for (final int afField : adFields) {
            totalReadCount += afField;
        }

        return new AlleleFrequencyData(alleleReadCount, totalReadCount);
    }
}
