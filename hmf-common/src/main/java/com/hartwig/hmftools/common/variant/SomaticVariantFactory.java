package com.hartwig.hmftools.common.variant;

import static com.hartwig.hmftools.common.variant.enrich.HighConfidenceEnrichment.HIGH_CONFIDENCE_FLAG;
import static com.hartwig.hmftools.common.variant.enrich.KataegisEnrichment.KATAEGIS_FLAG;
import static com.hartwig.hmftools.common.variant.enrich.PurityEnrichment.PURPLE_AF_INFO;
import static com.hartwig.hmftools.common.variant.enrich.PurityEnrichment.PURPLE_BIALLELIC_FLAG;
import static com.hartwig.hmftools.common.variant.enrich.PurityEnrichment.PURPLE_CN_INFO;
import static com.hartwig.hmftools.common.variant.enrich.PurityEnrichment.PURPLE_GERMLINE_INFO;
import static com.hartwig.hmftools.common.variant.enrich.PurityEnrichment.PURPLE_MINOR_ALLELE_PLOIDY_INFO;
import static com.hartwig.hmftools.common.variant.enrich.PurityEnrichment.PURPLE_PLOIDY_INFO;
import static com.hartwig.hmftools.common.variant.enrich.SomaticRefContextEnrichment.MICROHOMOLOGY_FLAG;
import static com.hartwig.hmftools.common.variant.enrich.SomaticRefContextEnrichment.REPEAT_COUNT_FLAG;
import static com.hartwig.hmftools.common.variant.enrich.SomaticRefContextEnrichment.REPEAT_SEQUENCE_FLAG;
import static com.hartwig.hmftools.common.variant.enrich.SomaticRefContextEnrichment.TRINUCLEOTIDE_FLAG;

import static htsjdk.tribble.AbstractFeatureReader.getFeatureReader;

import java.io.IOException;
import java.util.Arrays;
import java.util.List;
import java.util.Optional;
import java.util.StringJoiner;
import java.util.stream.Collectors;

import com.google.common.annotations.VisibleForTesting;
import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.purple.region.GermlineStatus;
import com.hartwig.hmftools.common.variant.cosmic.CosmicAnnotation;
import com.hartwig.hmftools.common.variant.cosmic.CosmicAnnotationFactory;
import com.hartwig.hmftools.common.variant.enrich.HotspotEnrichment;
import com.hartwig.hmftools.common.variant.enrich.NoSomaticEnrichment;
import com.hartwig.hmftools.common.variant.enrich.SomaticEnrichment;
import com.hartwig.hmftools.common.variant.enrich.SubclonalLikelihoodEnrichment;
import com.hartwig.hmftools.common.variant.enrich.VariantContextEnrichment;
import com.hartwig.hmftools.common.variant.enrich.VariantContextEnrichmentFactory;
import com.hartwig.hmftools.common.variant.filter.AlwaysPassFilter;
import com.hartwig.hmftools.common.variant.filter.ChromosomeFilter;
import com.hartwig.hmftools.common.variant.filter.NTFilter;
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

        return new SomaticVariantFactory(filter, enrichment, VariantContextEnrichmentFactory.noEnrichment());
    }

    @NotNull
    public static SomaticVariantFactory passOnlyInstance() {
        return filteredInstance(new PassingVariantFilter());
    }

    @NotNull
    public static SomaticVariantFactory passOnlyInstance(@NotNull VariantContextEnrichmentFactory factory) {
        return new SomaticVariantFactory(new PassingVariantFilter(), new NoSomaticEnrichment(), factory);
    }

    @NotNull
    public static SomaticVariantFactory filteredInstance(@NotNull VariantContextFilter... filters) {
        final CompoundFilter filter = new CompoundFilter(true);
        filter.addAll(Arrays.asList(filters));
        final SomaticEnrichment noEnrichment = new NoSomaticEnrichment();

        return new SomaticVariantFactory(filter, noEnrichment, VariantContextEnrichmentFactory.noEnrichment());
    }

    @NotNull
    public static SomaticVariantFactory filteredInstanceWithEnrichment(@NotNull VariantContextFilter filter,
            @NotNull SomaticEnrichment somaticEnrichment, @NotNull VariantContextEnrichmentFactory factory) {
        return new SomaticVariantFactory(filter, somaticEnrichment, factory);
    }

    private static final String DBSNP_IDENTIFIER = "rs";
    private static final String COSMIC_IDENTIFIER = "COSM";
    private static final String ID_SEPARATOR = ";";
    private static final String MAPPABILITY_TAG = "MAPPABILITY";
    private static final String RECOVERED_FLAG = "RECOVERED";

    public static final String PASS_FILTER = "PASS";
    private static final String NEAR_INDEL_PON_FILTER = "NEAR_INDEL_PON";

    @NotNull
    private final CompoundFilter filter;
    @NotNull
    private final SomaticEnrichment enrichment;
    @NotNull
    private final CanonicalAnnotation canonicalAnnotationFactory;
    @NotNull
    private final VariantContextEnrichmentFactory variantContextEnrichmentFactory;

    private SomaticVariantFactory(@NotNull final VariantContextFilter filter, @NotNull final SomaticEnrichment enrichment,
            @NotNull final VariantContextEnrichmentFactory variantContextEnrichmentFactory) {
        this.filter = new CompoundFilter(true);
        this.filter.add(new ChromosomeFilter());
        this.filter.add(new NTFilter());
        this.filter.add(filter);

        this.enrichment = enrichment;
        this.canonicalAnnotationFactory = new CanonicalAnnotation();
        this.variantContextEnrichmentFactory = variantContextEnrichmentFactory;

    }

    @NotNull
    public List<SomaticVariant> fromVCFFile(@NotNull final String sample, @NotNull final String vcfFile) throws IOException {
        final List<VariantContext> variants = Lists.newArrayList();
        final VariantContextEnrichment enrichment = variantContextEnrichmentFactory.create(variants::add);

        try (final AbstractFeatureReader<VariantContext, LineIterator> reader = getFeatureReader(vcfFile, new VCFCodec(), false)) {
            final VCFHeader header = (VCFHeader) reader.getHeader();
            if (!sampleInFile(sample, header)) {
                throw new IllegalArgumentException("Sample " + sample + " not found in vcf file " + vcfFile);
            }

            if (!header.hasFormatLine("AD")) {
                throw new IllegalArgumentException("Allelic depths is a required format field in vcf file " + vcfFile);
            }

            for (VariantContext variant : reader.iterator()) {
                // Note we need pon filtered indels for near indel pon logic to work correctly
                if (filter.test(variant) || NearPonFilteredIndel.isPonFilteredIndel(variant)) {
                    enrichment.accept(variant);
                }
            }

            enrichment.flush();
        }

        return process(sample, variants);
    }

    @NotNull
    private List<SomaticVariant> process(@NotNull final String sample, @NotNull final List<VariantContext> allVariantContexts) {
        final List<SomaticVariant> variants = Lists.newArrayList();

        for (int i = 0; i < allVariantContexts.size(); i++) {
            final VariantContext context = allVariantContexts.get(i);
            if (NearPonFilteredIndel.isNearPonFilteredIndel(i, allVariantContexts)) {
                context.getCommonInfo().addFilter(NEAR_INDEL_PON_FILTER);
            }

            createVariant(sample, context).ifPresent(variants::add);
        }

        return variants;
    }

    @NotNull
    @VisibleForTesting
    public Optional<SomaticVariant> createVariant(@NotNull final String sample, @NotNull final VariantContext context) {
        final Genotype genotype = context.getGenotype(sample);

        if (filter.test(context) && genotype.hasAD() && genotype.getAD().length > 1) {
            final AllelicDepth allelicDepth = AllelicDepth.fromGenotype(context.getGenotype(sample));
            if (allelicDepth.totalReadCount() > 0) {
                return Optional.of(createVariantBuilder(allelicDepth, context, canonicalAnnotationFactory))
                        .map(x -> enrichment.enrich(x, context))
                        .map(ImmutableSomaticVariantImpl.Builder::build);
            }
        }
        return Optional.empty();
    }

    @NotNull
    public SomaticVariant createSomaticVariant(@NotNull final String sample, @NotNull final VariantContext context) {
        final AllelicDepth allelicDepth = AllelicDepth.fromGenotype(context.getGenotype(sample));

        return Optional.of(createVariantBuilder(allelicDepth, context, canonicalAnnotationFactory))
                .map(x -> enrichment.enrich(x, context))
                .map(ImmutableSomaticVariantImpl.Builder::build)
                .get();
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
                .hotspot(HotspotEnrichment.fromVariant(context))
                .minorAllelePloidy(context.getAttributeAsDouble(PURPLE_MINOR_ALLELE_PLOIDY_INFO, 0))
                .adjustedCopyNumber(context.getAttributeAsDouble(PURPLE_CN_INFO, 0))
                .adjustedVAF(context.getAttributeAsDouble(PURPLE_AF_INFO, 0))
                .germlineStatus(GermlineStatus.valueOf(context.getAttributeAsString(PURPLE_GERMLINE_INFO, "UNKNOWN")))
                .ploidy(context.getAttributeAsDouble(PURPLE_PLOIDY_INFO, 0))
                .mappability(context.getAttributeAsDouble(MAPPABILITY_TAG, 0))
                .kataegis(context.getAttributeAsString(KATAEGIS_FLAG, Strings.EMPTY))
                .trinucleotideContext(context.getAttributeAsString(TRINUCLEOTIDE_FLAG, Strings.EMPTY))
                .microhomology(context.getAttributeAsString(MICROHOMOLOGY_FLAG, Strings.EMPTY))
                .repeatCount(context.getAttributeAsInt(REPEAT_COUNT_FLAG, 0))
                .repeatSequence(context.getAttributeAsString(REPEAT_SEQUENCE_FLAG, Strings.EMPTY))
                .subclonalLikelihood(context.getAttributeAsDouble(SubclonalLikelihoodEnrichment.SUBCLONAL_LIKELIHOOD_FLAG, 0))
                // Note: getAttributeAsBoolean(x, false) is safer than hasAttribute(x)
                .recovered(context.getAttributeAsBoolean(RECOVERED_FLAG, false))
                .biallelic(context.getAttributeAsBoolean(PURPLE_BIALLELIC_FLAG, false))
                .highConfidenceRegion(context.getAttributeAsBoolean(HIGH_CONFIDENCE_FLAG, false));

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

}
