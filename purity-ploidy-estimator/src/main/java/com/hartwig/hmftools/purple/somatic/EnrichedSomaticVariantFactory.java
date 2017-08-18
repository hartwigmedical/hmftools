package com.hartwig.hmftools.purple.somatic;

import static com.hartwig.hmftools.common.variant.ImmutableEnrichedSomaticVariant.Builder;
import static com.hartwig.hmftools.common.variant.ImmutableEnrichedSomaticVariant.builder;

import java.util.Collections;
import java.util.List;
import java.util.Optional;
import java.util.stream.Collectors;

import com.google.common.collect.Multimap;
import com.hartwig.hmftools.common.purple.PurityAdjuster;
import com.hartwig.hmftools.common.purple.copynumber.PurpleCopyNumber;
import com.hartwig.hmftools.common.purple.gender.Gender;
import com.hartwig.hmftools.common.purple.purity.FittedPurity;
import com.hartwig.hmftools.common.purple.repeat.RepeatContextFactory;
import com.hartwig.hmftools.common.region.GenomeRegion;
import com.hartwig.hmftools.common.region.GenomeRegionSelector;
import com.hartwig.hmftools.common.region.GenomeRegionSelectorFactory;
import com.hartwig.hmftools.common.variant.EnrichedSomaticVariant;
import com.hartwig.hmftools.common.variant.SomaticVariant;
import com.hartwig.hmftools.common.variant.VariantAnnotation;

import org.jetbrains.annotations.NotNull;

import htsjdk.samtools.reference.IndexedFastaSequenceFile;
import htsjdk.samtools.reference.ReferenceSequence;

public class EnrichedSomaticVariantFactory {

    @NotNull
    private final PurityAdjuster purityAdjuster;
    @NotNull
    private final GenomeRegionSelector<GenomeRegion> highConfidenceSelector;
    @NotNull
    private final GenomeRegionSelector<PurpleCopyNumber> copyNumberSelector;
    @NotNull
    private final Optional<IndexedFastaSequenceFile> optionalReference;

    public EnrichedSomaticVariantFactory(double purity, double normFactor,
            @NotNull  final Multimap<String, GenomeRegion> highConfidenceRegions,
            @NotNull  final Multimap<String, PurpleCopyNumber> copyNumbers, @NotNull final IndexedFastaSequenceFile reference) {
        purityAdjuster = new PurityAdjuster(Gender.MALE, purity, normFactor);
        highConfidenceSelector = GenomeRegionSelectorFactory.create(highConfidenceRegions);
        copyNumberSelector = GenomeRegionSelectorFactory.create(copyNumbers);
        this.optionalReference = Optional.of(reference);
    }

    public EnrichedSomaticVariantFactory(@NotNull  final FittedPurity purity, @NotNull final List<PurpleCopyNumber> copyNumbers) {
        purityAdjuster = new PurityAdjuster(Gender.MALE, purity);
        highConfidenceSelector = GenomeRegionSelectorFactory.create(Collections.emptyList());
        copyNumberSelector = GenomeRegionSelectorFactory.create(copyNumbers);
        this.optionalReference = Optional.empty();
    }

    public List<EnrichedSomaticVariant> enrich(final List<SomaticVariant> variants) {
        return variants.stream().map(this::enrich).collect(Collectors.toList());
    }

    private EnrichedSomaticVariant enrich(@NotNull final SomaticVariant variant) {
        final Builder builder = createBuilder(variant);

        highConfidenceSelector.select(variant).ifPresent(x -> inHighConfidenceRegion(builder));
        copyNumberSelector.select(variant).ifPresent(x -> addCopyNumber(builder, x, variant.alleleFrequency()));
        addAnnotations(builder, variant);
        addTrinucleotideContext(builder, variant);
        addGenomeContext(builder, variant);

        return builder.build();
    }

    private Builder addAnnotations(@NotNull final Builder builder, @NotNull final SomaticVariant variant) {
        final List<VariantAnnotation> annotations = variant.annotations();
        if (!annotations.isEmpty()) {
            VariantAnnotation annotation = annotations.get(0);
            builder.gene(annotation.gene());
            builder.effect(annotation.consequenceString());
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
                .adjustedVAF(0);
    }

    private Builder addGenomeContext(@NotNull final Builder builder, @NotNull final SomaticVariant variant) {
        if (optionalReference.isPresent()) {
            final IndexedFastaSequenceFile reference = optionalReference.get();
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
        }
        return builder;
    }

    private Builder addCopyNumber(@NotNull final Builder builder, @NotNull final PurpleCopyNumber copyNumber, double alleleFrequency) {
        double adjustedVAF = purityAdjuster.purityAdjustedVAF(Math.max(0.001, copyNumber.averageTumorCopyNumber()), alleleFrequency);
        return builder.adjustedCopyNumber(copyNumber.averageTumorCopyNumber()).adjustedVAF(adjustedVAF);
    }

    private Builder addTrinucleotideContext(@NotNull final Builder builder, @NotNull final SomaticVariant variant) {
        if (optionalReference.isPresent()) {
            final IndexedFastaSequenceFile reference = optionalReference.get();
            final ReferenceSequence sequence =
                    reference.getSubsequenceAt(variant.chromosome(), Math.max(1, variant.position() - 1), variant.position() + 1);
            return builder.trinucleotideContext(sequence.getBaseString());
        }

        return builder;
    }

    private Builder inHighConfidenceRegion(@NotNull final Builder builder) {
        return builder.highConfidenceRegion(true);
    }
}
