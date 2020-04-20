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
import com.hartwig.hmftools.common.variant.snpeff.SnpEffSummary;
import com.hartwig.hmftools.common.variant.snpeff.SnpEffSummaryFactory;

import org.apache.logging.log4j.util.Strings;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

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

    static final String PASS_FILTER = "PASS";
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

    public List<SomaticVariant> fromVCFFile(@NotNull final String tumor, @NotNull final String vcfFile) throws IOException {
        return fromVCFFile(tumor, null, null, vcfFile);
    }

    @NotNull
    public List<SomaticVariant> fromVCFFile(@NotNull final String tumor, @Nullable final String reference, @Nullable final String rna,
            @NotNull final String vcfFile) throws IOException {
        final List<VariantContext> variants = Lists.newArrayList();
        final VariantContextEnrichment enrichment = variantContextEnrichmentFactory.create(variants::add);

        try (final AbstractFeatureReader<VariantContext, LineIterator> reader = getFeatureReader(vcfFile, new VCFCodec(), false)) {
            final VCFHeader header = (VCFHeader) reader.getHeader();
            if (!sampleInFile(tumor, header)) {
                throw new IllegalArgumentException("Sample " + tumor + " not found in vcf file " + vcfFile);
            }

            if (reference != null && !sampleInFile(reference, header)) {
                throw new IllegalArgumentException("Sample " + reference + " not found in vcf file " + vcfFile);
            }

            if (rna != null && !sampleInFile(rna, header)) {
                throw new IllegalArgumentException("Sample " + rna + " not found in vcf file " + vcfFile);
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

        return process(tumor, reference, rna, variants);
    }

    @NotNull
    private List<SomaticVariant> process(@NotNull final String sample, @Nullable final String reference, @Nullable final String rna,
            @NotNull final List<VariantContext> allVariantContexts) {
        final List<SomaticVariant> variants = Lists.newArrayList();

        for (int i = 0; i < allVariantContexts.size(); i++) {
            final VariantContext context = allVariantContexts.get(i);
            if (NearPonFilteredIndel.isNearPonFilteredIndel(i, allVariantContexts)) {
                context.getCommonInfo().addFilter(NEAR_INDEL_PON_FILTER);
            }

            createVariant(sample, reference, rna, context).ifPresent(variants::add);
        }

        return variants;
    }

    @NotNull
    public Optional<SomaticVariant> createVariant(@NotNull final String sample, @NotNull final VariantContext context) {
        return createVariant(sample, null, null, context);
    }

    @NotNull
    public Optional<SomaticVariant> createVariant(@NotNull final String sample, @Nullable final String reference,
            @Nullable final String rna, @NotNull final VariantContext context) {
        final Genotype genotype = context.getGenotype(sample);

        if (filter.test(context) && AllelicDepth.containsAllelicDepth(genotype)) {
            final AllelicDepth tumorDepth = AllelicDepth.fromGenotype(context.getGenotype(sample));

            final Optional<AllelicDepth> referenceDepth = Optional.ofNullable(reference)
                    .flatMap(x -> Optional.ofNullable(context.getGenotype(x)))
                    .filter(AllelicDepth::containsAllelicDepth)
                    .map(AllelicDepth::fromGenotype);

            final Optional<AllelicDepth> rnaDepth = Optional.ofNullable(rna)
                    .flatMap(x -> Optional.ofNullable(context.getGenotype(x)))
                    .filter(AllelicDepth::containsAllelicDepth)
                    .map(AllelicDepth::fromGenotype);

            if (tumorDepth.totalReadCount() > 0) {
                return Optional.of(createVariantBuilder(tumorDepth, context, canonicalAnnotationFactory))
                        .map(x -> x.rnaDepth(rnaDepth.orElse(null)))
                        .map(x -> x.referenceDepth(referenceDepth.orElse(null)))
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
        SnpEffSummaryFactory snpEffSummaryFactory = new SnpEffSummaryFactory(canonicalAnnotationFactory);

        ImmutableSomaticVariantImpl.Builder builder = ImmutableSomaticVariantImpl.builder()
                .qual(context.getPhredScaledQual())
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
                .tier(VariantTier.fromString(context.getAttributeAsString("TIER", VariantTier.UNKNOWN.toString())))
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
        attachSnpEffAnnotations(builder, context, snpEffSummaryFactory);
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
        final Optional<CosmicAnnotation> canonicalCosmicAnnotation = canonicalAnnotationFactory.canonicalCosmicAnnotation(cosmicAnnotations);

        if (canonicalCosmicAnnotation.isPresent()) {
            builder.canonicalCosmicID(canonicalCosmicAnnotation.get().id());
        } else if (!cosmicIDs.isEmpty()) {
            builder.canonicalCosmicID(cosmicIDs.get(0));
        }
    }

    private static void attachSnpEffAnnotations(@NotNull final ImmutableSomaticVariantImpl.Builder builder, @NotNull VariantContext context,
            @NotNull SnpEffSummaryFactory snpEffSummaryFactory) {
        final SnpEffSummary snpEffSummary = snpEffSummaryFactory.fromAnnotations(context);

        builder.worstEffect(snpEffSummary.worstEffect())
                .worstCodingEffect(snpEffSummary.worstCodingEffect())
                .worstEffectTranscript(snpEffSummary.worstTranscript())
                .canonicalEffect(snpEffSummary.canonicalEffect())
                .canonicalCodingEffect(snpEffSummary.canonicalCodingEffect())
                .canonicalHgvsCodingImpact(snpEffSummary.canonicalHgvsCodingImpact())
                .canonicalHgvsProteinImpact(snpEffSummary.canonicalHgvsProteinImpact())
                .gene(snpEffSummary.gene())
                .genesAffected(snpEffSummary.genesAffected());
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
