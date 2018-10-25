package com.hartwig.hmftools.common.variant;

import static htsjdk.tribble.AbstractFeatureReader.getFeatureReader;

import java.io.IOException;
import java.util.Arrays;
import java.util.List;
import java.util.Optional;
import java.util.StringJoiner;
import java.util.stream.Collectors;

import com.google.common.annotations.VisibleForTesting;
import com.google.common.base.Preconditions;
import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.variant.cosmic.CosmicAnnotation;
import com.hartwig.hmftools.common.variant.cosmic.CosmicAnnotationFactory;
import com.hartwig.hmftools.common.variant.enrich.NoSomaticEnrichment;
import com.hartwig.hmftools.common.variant.enrich.SomaticEnrichment;
import com.hartwig.hmftools.common.variant.filter.AlwaysPassFilter;
import com.hartwig.hmftools.common.variant.filter.ChromosomeFilter;
import com.hartwig.hmftools.common.variant.filter.HotspotFilter;
import com.hartwig.hmftools.common.variant.filter.NearIndelPonFilter;
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

    @NotNull
    public static SomaticVariantFactory unfilteredInstance() {
        final VariantContextFilter filter = new AlwaysPassFilter();
        final SomaticEnrichment enrichment = new NoSomaticEnrichment();

        return new SomaticVariantFactory(filter, enrichment);
    }

    @NotNull
    public static SomaticVariantFactory passOnlyInstance() {
        return filteredInstance(new PassingVariantFilter());
    }

    @NotNull
    public static SomaticVariantFactory filteredInstance(@NotNull VariantContextFilter... filters) {
        final CompoundFilter filter = new CompoundFilter(true);
        filter.addAll(Arrays.asList(filters));
        final SomaticEnrichment noEnrichment = new NoSomaticEnrichment();

        return new SomaticVariantFactory(filter, noEnrichment);
    }

    @NotNull
    public static SomaticVariantFactory filteredInstanceWithEnrichment(@NotNull VariantContextFilter filter,
            @NotNull SomaticEnrichment somaticEnrichment) {
        return new SomaticVariantFactory(filter, somaticEnrichment);
    }

    private static final HotspotFilter HOTSPOT_FILTER = new HotspotFilter();

    private static final String DBSNP_IDENTIFIER = "rs";
    private static final String COSMIC_IDENTIFIER = "COSM";
    private static final String ID_SEPARATOR = ";";
    private static final String MAPPABILITY_TAG = "MAPPABILITY";

    static final String PASS_FILTER = "PASS";
    private static final String NEAR_INDEL_PON_FILTER = "NEAR_INDEL_PON";

    @NotNull
    private final CompoundFilter filter;
    @NotNull
    private final SomaticEnrichment enrichment;
    @NotNull
    private final CanonicalAnnotation canonicalAnnotationFactory;

    private SomaticVariantFactory(@NotNull final VariantContextFilter filter, @NotNull final SomaticEnrichment enrichment) {
        this.filter = new CompoundFilter(true);
        this.filter.add(new ChromosomeFilter());
        this.filter.add(filter);

        this.enrichment = enrichment;
        this.canonicalAnnotationFactory = new CanonicalAnnotation();
    }

    @NotNull
    public List<SomaticVariant> fromVCFFile(@NotNull final String sample, @NotNull final String vcfFile) throws IOException {
        final List<VariantContext> variants = Lists.newArrayList();

        try (final AbstractFeatureReader<VariantContext, LineIterator> reader = getFeatureReader(vcfFile, new VCFCodec(), false)) {
            final VCFHeader header = (VCFHeader) reader.getHeader();
            if (!sampleInFile(sample, header)) {
                throw new IllegalArgumentException("Sample " + sample + " not found in vcf file " + vcfFile);
            }

            if (!header.hasFormatLine("AD")) {
                throw new IllegalArgumentException("Allelic depths is a required format field in vcf file " + vcfFile);
            }

            variants.addAll(reader.iterator().toList());
        }

        return process(sample, variants);
    }

    @VisibleForTesting
    List<SomaticVariant> process(@NotNull final String sample, @NotNull final List<VariantContext> allVariantContexts) {
        final List<SomaticVariant> variants = Lists.newArrayList();

        for (int i = 0; i < allVariantContexts.size(); i++) {
            final VariantContext context = allVariantContexts.get(i);
            if (NearIndelPonFilter.isIndelNearPon(i, allVariantContexts)) {
                context.getCommonInfo().addFilter(NEAR_INDEL_PON_FILTER);
            }

            createVariant(sample, context).ifPresent(variants::add);
        }

        return variants;
    }

    @NotNull
    // TODO (KODU): This function is used by BachelorPP, should probably change.
    public Optional<SomaticVariant> createVariant(@NotNull final String sample, @NotNull final VariantContext context) {
        final Genotype genotype = context.getGenotype(sample);

        if (filter.test(context) && genotype.hasAD() && genotype.getAD().length > 1) {
            final AllelicDepth allelicDepth = determineAlleleFrequencies(context.getGenotype(sample));
            if (allelicDepth.totalReadCount() > 0) {
                return Optional.of(createVariantBuilder(allelicDepth, context, canonicalAnnotationFactory))
                        .map(x -> enrichment.enrich(x, context))
                        .map(ImmutableSomaticVariantImpl.Builder::build);
            }
        }
        return Optional.empty();
    }

    @NotNull
    private static ImmutableSomaticVariantImpl.Builder createVariantBuilder(@NotNull final AllelicDepth allelicDepth,
            @NotNull final VariantContext context, @NotNull CanonicalAnnotation canonicalAnnotationFactory) {

        ImmutableSomaticVariantImpl.Builder builder = ImmutableSomaticVariantImpl.builder()
                .chromosome(context.getContig())
                .position(context.getStart())
                .ref(context.getReference().getBaseString())
                .alt(alt(context))
                .alleleReadCount(allelicDepth.alleleReadCount())
                .totalReadCount(allelicDepth.totalReadCount())
                .hotspot(HOTSPOT_FILTER.test(context) ? Hotspot.HOTSPOT : Hotspot.NON_HOTSPOT)
                .mappability(context.getAttributeAsDouble(MAPPABILITY_TAG, 0));

        attachIDAndCosmicAnnotations(builder, context, canonicalAnnotationFactory);
        attachSnpEffAnnotations(builder, context, canonicalAnnotationFactory);
        attachFilter(builder, context);
        attachType(builder, context);

        return builder;
    }

    private static void attachIDAndCosmicAnnotations(@NotNull final ImmutableSomaticVariantImpl.Builder builder,
            @NotNull VariantContext context, @NotNull CanonicalAnnotation canonicalAnnotationFactory) {
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
                canonicalAnnotationFactory.canonicalCosmicAnnotation(cosmicAnnotations);

        if (canonicalCosmicAnnotation.isPresent()) {
            builder.canonicalCosmicID(canonicalCosmicAnnotation.get().id());
        } else if (!cosmicIDs.isEmpty()) {
            builder.canonicalCosmicID(cosmicIDs.get(0));
        }
    }

    private static void attachSnpEffAnnotations(@NotNull final ImmutableSomaticVariantImpl.Builder builder, @NotNull VariantContext context,
            @NotNull CanonicalAnnotation canonicalAnnotationFactory) {
        final List<SnpEffAnnotation> allAnnotations = SnpEffAnnotationFactory.fromContext(context);
        builder.snpEffAnnotations(allAnnotations);

        final List<SnpEffAnnotation> transcriptAnnotations =
                allAnnotations.stream().filter(SnpEffAnnotation::isTranscriptFeature).collect(Collectors.toList());
        if (!transcriptAnnotations.isEmpty()) {
            final SnpEffAnnotation worstAnnotation = transcriptAnnotations.get(0);
            builder.worstEffect(worstAnnotation.consequenceString());
            builder.worstCodingEffect(CodingEffect.effect(worstAnnotation.gene(), worstAnnotation.consequences()));
            builder.worstEffectTranscript(worstAnnotation.transcript());
        } else {
            builder.worstEffect(Strings.EMPTY);
            builder.worstCodingEffect(CodingEffect.UNDEFINED);
            builder.worstEffectTranscript(Strings.EMPTY);
        }

        final Optional<SnpEffAnnotation> canonicalAnnotation = canonicalAnnotationFactory.canonicalSnpEffAnnotation(transcriptAnnotations);
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
        final int alleleReadCount = adFields[1];
        int totalReadCount = 0;
        for (final int afField : adFields) {
            totalReadCount += afField;
        }

        return ImmutableAllelicDepthImpl.builder().alleleReadCount(alleleReadCount).totalReadCount(totalReadCount).build();
    }
}
