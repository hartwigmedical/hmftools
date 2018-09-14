package com.hartwig.hmftools.common.variant;

import static htsjdk.tribble.AbstractFeatureReader.getFeatureReader;

import java.io.IOException;
import java.util.Arrays;
import java.util.List;
import java.util.Map;
import java.util.Optional;
import java.util.StringJoiner;
import java.util.stream.Collectors;

import com.google.common.annotations.VisibleForTesting;
import com.google.common.base.Preconditions;
import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.gene.CanonicalTranscriptFactory;
import com.hartwig.hmftools.common.gene.TranscriptRegion;
import com.hartwig.hmftools.common.variant.cosmic.CosmicAnnotation;
import com.hartwig.hmftools.common.variant.cosmic.CosmicAnnotationFactory;
import com.hartwig.hmftools.common.variant.filter.ChromosomeFilter;
import com.hartwig.hmftools.common.variant.filter.HotspotFilter;
import com.hartwig.hmftools.common.variant.snpeff.SnpEffAnnotation;
import com.hartwig.hmftools.common.variant.snpeff.SnpEffAnnotationFactory;

import org.apache.logging.log4j.util.Strings;
import org.jetbrains.annotations.NotNull;

import htsjdk.tribble.AbstractFeatureReader;
import htsjdk.tribble.readers.LineIterator;
import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.Genotype;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.filter.CompoundFilter;
import htsjdk.variant.variantcontext.filter.PassingVariantFilter;
import htsjdk.variant.variantcontext.filter.VariantContextFilter;
import htsjdk.variant.vcf.VCFCodec;
import htsjdk.variant.vcf.VCFHeader;

public class SomaticVariantFactory {

    static final String PASS_FILTER = "PASS";

    @NotNull
    public static SomaticVariantFactory unfilteredInstance() {
        return filteredInstance();
    }

    @NotNull
    public static SomaticVariantFactory passOnlyInstance() {
        return filteredInstance(new PassingVariantFilter());
    }

    @NotNull
    public static SomaticVariantFactory filteredInstance(@NotNull VariantContextFilter... filters) {
        final CompoundFilter filter = new CompoundFilter(true);
        filter.add(new ChromosomeFilter());
        filter.addAll(Arrays.asList(filters));
        return new SomaticVariantFactory(filter);
    }

    private static final HotspotFilter HOTSPOT_FILTER = new HotspotFilter();

    private static final String DBSNP_IDENTIFIER = "rs";
    private static final String COSMIC_IDENTIFIER = "COSM";
    private static final String ID_SEPARATOR = ";";
    private static final String MAPPABILITY_TAG = "MAPPABILITY";

    @NotNull
    private final VariantContextFilter filter;
    @NotNull
    private final Map<String, String> transcriptIdToGeneMap;

    private SomaticVariantFactory(@NotNull final VariantContextFilter filter) {
        this.filter = filter;
        this.transcriptIdToGeneMap = CanonicalTranscriptFactory.create()
                .stream()
                .collect(Collectors.toMap(TranscriptRegion::transcriptID, TranscriptRegion::gene));
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
    public Optional<SomaticVariant> createVariant(@NotNull final String sample, @NotNull final VariantContext context) {
        if (filter.test(context)) {
            final Genotype genotype = context.getGenotype(sample);
            if (genotype.hasAD() && genotype.getAD().length > 1) {
                final AllelicDepth frequencyData = determineAlleleFrequencies(genotype);
                if (frequencyData.totalReadCount() > 0) {
                    ImmutableSomaticVariantImpl.Builder builder = ImmutableSomaticVariantImpl.builder()
                            .chromosome(context.getContig())
                            .position(context.getStart())
                            .ref(context.getReference().getBaseString())
                            .alt(alt(context))
                            .alleleReadCount(frequencyData.alleleReadCount())
                            .totalReadCount(frequencyData.totalReadCount())
                            .hotspot(HOTSPOT_FILTER.test(context))
                            .mappability(context.getAttributeAsDouble(MAPPABILITY_TAG, 0));

                    attachIDAndCosmicAnnotations(builder, context);
                    attachSnpEffAnnotations(builder, context);
                    attachFilter(builder, context);
                    attachType(builder, context);
                    return Optional.of(builder.build());
                }
            }
        }
        return Optional.empty();
    }

    private void attachIDAndCosmicAnnotations(@NotNull final ImmutableSomaticVariantImpl.Builder builder, @NotNull VariantContext context) {
        final String ID = context.getID();
        final List<String> cosmicIDs = Lists.newArrayList();
        if (!ID.isEmpty()) {
            final String[] ids = ID.split(ID_SEPARATOR);
            for (final String id : ids) {
                if (id.contains(DBSNP_IDENTIFIER)) {
                    builder.dbsnpID(id);
                } else if (id.contains(COSMIC_IDENTIFIER)) {
                    cosmicIDs.add(id);
                }
            }
        }
        builder.cosmicIDs(cosmicIDs);

        final List<CosmicAnnotation> cosmicAnnotations = CosmicAnnotationFactory.fromContext(context);
        builder.cosmicAnnotations(cosmicAnnotations);

        final Optional<CosmicAnnotation> canonicalCosmicAnnotation =
                cosmicAnnotations.stream().filter(x -> transcriptIdToGeneMap.keySet().contains(x.transcript())).findFirst();

        if (canonicalCosmicAnnotation.isPresent()) {
            builder.canonicalCosmicID(canonicalCosmicAnnotation.get().id());
        } else if (!cosmicIDs.isEmpty()) {
            builder.canonicalCosmicID(cosmicIDs.get(0));
        }
    }

    private void attachSnpEffAnnotations(@NotNull final ImmutableSomaticVariantImpl.Builder builder, @NotNull VariantContext context) {
        final List<SnpEffAnnotation> allAnnotations = SnpEffAnnotationFactory.fromContext(context);
        builder.snpEffAnnotations(allAnnotations);

        final List<SnpEffAnnotation> transcriptAnnotations =
                allAnnotations.stream().filter(SnpEffAnnotation::isTranscriptFeature).collect(Collectors.toList());
        if (!transcriptAnnotations.isEmpty()) {
            final SnpEffAnnotation snpEffAnnotation = transcriptAnnotations.get(0);
            builder.worstEffect(snpEffAnnotation.consequenceString());
            builder.worstCodingEffect(CodingEffect.effect(snpEffAnnotation.gene(), snpEffAnnotation.consequences()));
            builder.worstEffectTranscript(snpEffAnnotation.transcript());
        } else {
            builder.worstEffect(Strings.EMPTY);
            builder.worstCodingEffect(CodingEffect.UNDEFINED);
            builder.worstEffectTranscript(Strings.EMPTY);
        }

        final Optional<SnpEffAnnotation> canonicalAnnotation =
                transcriptAnnotations.stream().filter(x -> transcriptIdToGeneMap.keySet().contains(x.transcript())).findFirst();
        if (canonicalAnnotation.isPresent()) {
            final SnpEffAnnotation annotation = canonicalAnnotation.get();
            builder.canonicalEffect(annotation.consequenceString());
            builder.canonicalCodingEffect(CodingEffect.effect(annotation.gene(), annotation.consequences()));
            builder.canonicalHgvsCodingImpact(annotation.hgvsCoding());
            builder.canonicalHgvsProteinImpact(annotation.hgvsProtein());
        } else {
            builder.canonicalEffect(Strings.EMPTY);
            builder.canonicalCodingEffect(CodingEffect.UNDEFINED);
            builder.canonicalHgvsCodingImpact(Strings.EMPTY);
            builder.canonicalHgvsProteinImpact(Strings.EMPTY);
        }

        final String firstGene = transcriptAnnotations.isEmpty() ? Strings.EMPTY : transcriptAnnotations.get(0).gene();
        final String gene = canonicalAnnotation.map(SnpEffAnnotation::gene).orElse(firstGene);
        builder.gene(gene);

        builder.genesEffected((int) transcriptAnnotations.stream()
                .map(SnpEffAnnotation::gene)
                .filter(x -> !x.isEmpty())
                .distinct()
                .count());
    }

    private static void attachFilter(@NotNull final ImmutableSomaticVariantImpl.Builder builder, @NotNull VariantContext context) {
        if (context.isFiltered()) {
            StringJoiner joiner = new StringJoiner(";");
            context.getFilters().forEach(joiner::add);
            builder.filter(joiner.toString());
        } else {
            builder.filter(PASS_FILTER);
        }
    }

    @NotNull
    private static VariantType type(@NotNull VariantContext context) {
        switch (context.getType()) {
            case MNP:
                return VariantType.MNP;
            case SNP:
                return VariantType.SNP;
            case INDEL:
                return VariantType.INDEL;
        }
        return VariantType.UNDEFINED;
    }

    private static void attachType(@NotNull final ImmutableSomaticVariantImpl.Builder builder, @NotNull VariantContext context) {
        builder.type(type(context));
    }

    private static boolean sampleInFile(@NotNull final String sample, @NotNull final VCFHeader header) {
        return header.getSampleNamesInOrder().stream().anyMatch(x -> x.equals(sample));
    }

    @NotNull
    private static String alt(@NotNull final VariantContext context) {
        return String.join(",", context.getAlternateAlleles().stream().map(Allele::toString).collect(Collectors.toList()));
    }

    @NotNull
    private static AllelicDepth determineAlleleFrequencies(@NotNull final Genotype genotype) {
        Preconditions.checkArgument(genotype.hasAD());

        int[] adFields = genotype.getAD();
        int totalReadCount = 0;
        final int alleleReadCount = adFields[1];
        for (final int afField : adFields) {
            totalReadCount += afField;
        }

        return ImmutableAllelicDepthImpl.builder().alleleReadCount(alleleReadCount).totalReadCount(totalReadCount).build();
    }
}
