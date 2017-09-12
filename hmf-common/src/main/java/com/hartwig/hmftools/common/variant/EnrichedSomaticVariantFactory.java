package com.hartwig.hmftools.common.variant;

import static com.hartwig.hmftools.common.variant.ImmutableEnrichedSomaticVariant.Builder;
import static com.hartwig.hmftools.common.variant.ImmutableEnrichedSomaticVariant.builder;

import java.util.List;
import java.util.Optional;
import java.util.stream.Collectors;

import com.google.common.collect.Multimap;
import com.hartwig.hmftools.common.gc.GCMedianReadCount;
import com.hartwig.hmftools.common.gc.GCProfile;
import com.hartwig.hmftools.common.purple.PurityAdjuster;
import com.hartwig.hmftools.common.purple.copynumber.PurpleCopyNumber;
import com.hartwig.hmftools.common.purple.gender.Gender;
import com.hartwig.hmftools.common.purple.repeat.RepeatContextFactory;
import com.hartwig.hmftools.common.region.GenomeRegion;
import com.hartwig.hmftools.common.region.GenomeRegionSelector;
import com.hartwig.hmftools.common.region.GenomeRegionSelectorFactory;

import org.apache.commons.math3.distribution.BinomialDistribution;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.jetbrains.annotations.NotNull;

import htsjdk.samtools.reference.IndexedFastaSequenceFile;
import htsjdk.samtools.reference.ReferenceSequence;

public class EnrichedSomaticVariantFactory {

    private static final Logger LOGGER = LogManager.getLogger(EnrichedSomaticVariantFactory.class);

    private static final int TRIALS = 100_000;

    @NotNull
    private final PurityAdjuster purityAdjuster;
    @NotNull
    private final GenomeRegionSelector<GenomeRegion> highConfidenceSelector;
    @NotNull
    private final GenomeRegionSelector<PurpleCopyNumber> copyNumberSelector;
    @NotNull
    private final GenomeRegionSelector<GCProfile> gcProfileSelector;
    @NotNull
    private final GCMedianReadCount medianReadCount;
    @NotNull
    private final IndexedFastaSequenceFile reference;

    private final double monoploidProbability;
    private final BinomialDistribution monoploidDistribution;
    private int unmatchedAnnotations;

    public EnrichedSomaticVariantFactory(double purity, double normFactor, int medianVariantTotalReadCount,
            @NotNull final Multimap<String, GenomeRegion> highConfidenceRegions,
            @NotNull final Multimap<String, PurpleCopyNumber> copyNumbers, @NotNull final Multimap<String, GCProfile> gcProfiles,
            @NotNull final GCMedianReadCount medianReadCount, @NotNull final IndexedFastaSequenceFile reference) {
        purityAdjuster = new PurityAdjuster(Gender.MALE, purity, normFactor);
        highConfidenceSelector = GenomeRegionSelectorFactory.create(highConfidenceRegions);
        copyNumberSelector = GenomeRegionSelectorFactory.create(copyNumbers);
        gcProfileSelector = GenomeRegionSelectorFactory.create(gcProfiles);
        this.medianReadCount = medianReadCount;
        this.reference = reference;

        this.monoploidProbability = medianVariantTotalReadCount * purity / (purity * purityAdjuster.impliedPloidy() + 2 * (1 - purity));
        this.monoploidDistribution = new BinomialDistribution(TRIALS, monoploidProbability / TRIALS);
    }

    public List<EnrichedSomaticVariant> enrich(final List<SomaticVariant> variants) {

        final List<EnrichedSomaticVariant> result = variants.stream().map(this::enrich).collect(Collectors.toList());
        if (unmatchedAnnotations > 0) {
            LOGGER.warn("There were {} unmatched annotated genes.", unmatchedAnnotations);
        }

        return result;
    }

    private EnrichedSomaticVariant enrich(@NotNull final SomaticVariant variant) {
        final Builder builder = createBuilder(variant);

        highConfidenceSelector.select(variant).ifPresent(x -> inHighConfidenceRegion(builder));

        final Optional<GCProfile> optionalGCProfile = gcProfileSelector.select(variant);
        final Optional<PurpleCopyNumber> optionalCopyNumber = copyNumberSelector.select(variant);
        if (optionalCopyNumber.isPresent()) {
            final PurpleCopyNumber copyNumber = optionalCopyNumber.get();
            builder.purityAdjustment(purityAdjuster, copyNumber, variant);

            BinomialDistribution inconsistentDistribution =
                    new BinomialDistribution(TRIALS, Math.max(0, copyNumber.averageTumorCopyNumber()) * monoploidProbability / TRIALS);

            double scalingFactor = optionalGCProfile.map(x -> 1d * medianReadCount.medianReadCount(x) / medianReadCount.medianReadCount()).orElse(1d);
            long adjustedAlleleCount = Math.round(variant.alleleReadCount() * scalingFactor);

            builder.clonality(adjustedAlleleCount, monoploidDistribution, inconsistentDistribution);
        }

        addAnnotations(builder, variant);
        addTrinucleotideContext(builder, variant);
        addGenomeContext(builder, variant);

        return builder.build();
    }

    private Builder addAnnotations(@NotNull final Builder builder, @NotNull final SomaticVariant variant) {
        final List<VariantAnnotation> annotations = variant.annotations();
        if (!annotations.isEmpty()) {
            // MIVO: get the first annotation for now, eventually we will want all
            final VariantAnnotation variantAnnotation = variant.annotations().get(0);
            variant.annotations().forEach(annotation -> {
                if (!annotation.gene().equals(variantAnnotation.gene())) {
                    unmatchedAnnotations++;
                    LOGGER.debug("Annotated gene (" + annotation.gene() + ") does not match gene expected from first annotation ( "
                            + variantAnnotation.gene() + ") for variant: " + variant);
                }
            });
            builder.gene(variantAnnotation.gene());
            builder.effect(variantAnnotation.consequenceString());
        }
        return builder;
    }

    private Builder createBuilder(@NotNull final SomaticVariant variant) {
        String cosmicId = variant.cosmicID();
        String dbsnpId = variant.dbsnpID();

        return builder().from(variant)
                .trinucleotideContext("")
                .microhomology("")
                .refGenomeContext("")
                .gene("")
                .cosmicId(cosmicId == null ? "" : cosmicId)
                .dbsnpId(dbsnpId == null ? "" : dbsnpId)
                .effect("")
                .repeatCount(0)
                .repeatSequence("")
                .totalReadCount(variant.totalReadCount())
                .alleleReadCount(variant.alleleReadCount())
                .highConfidenceRegion(false)
                .adjustedCopyNumber(0)
                .clonality(Clonality.UNKNOWN)
                .lossOfHeterozygosity(false)
                .adjustedVAF(0);
    }

    private Builder addGenomeContext(@NotNull final Builder builder, @NotNull final SomaticVariant variant) {
        long positionBeforeEvent = variant.position();
        long start = Math.max(positionBeforeEvent - 100, 1);
        long maxEnd = reference.getSequenceDictionary().getSequence(variant.chromosome()).getSequenceLength() - 1;
        long end = Math.min(positionBeforeEvent + 100, maxEnd);
        int relativePosition = (int) (positionBeforeEvent - start);
        final String sequence = reference.getSubsequenceAt(variant.chromosome(), start, end).getBaseString();
        builder.refGenomeContext(sequence);

        RepeatContextFactory.repeats(relativePosition, sequence, variant.ref(), variant.alt())
                .ifPresent(x -> builder.repeatSequence(x.sequence()).repeatCount(x.count()));

        if (variant.ref().length() != variant.alt().length()) {
            final String microhomology = Microhomology.microhomology(relativePosition, sequence, variant.ref(), variant.alt());
            return builder.microhomology(microhomology);
        }
        return builder;
    }

    private Builder addTrinucleotideContext(@NotNull final Builder builder, @NotNull final SomaticVariant variant) {
        final ReferenceSequence sequence =
                reference.getSubsequenceAt(variant.chromosome(), Math.max(1, variant.position() - 1), variant.position() + 1);
        return builder.trinucleotideContext(sequence.getBaseString());
    }

    private Builder inHighConfidenceRegion(@NotNull final Builder builder) {
        return builder.highConfidenceRegion(true);
    }
}
