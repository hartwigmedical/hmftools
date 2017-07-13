package com.hartwig.hmftools.purple.somatic;

import static com.hartwig.hmftools.common.variant.ImmutableEnrichedSomaticVariant.Builder;
import static com.hartwig.hmftools.common.variant.ImmutableEnrichedSomaticVariant.builder;

import java.util.List;
import java.util.stream.Collectors;

import com.google.common.collect.Multimap;
import com.hartwig.hmftools.common.purple.PurityAdjuster;
import com.hartwig.hmftools.common.purple.copynumber.PurpleCopyNumber;
import com.hartwig.hmftools.common.purple.gender.Gender;
import com.hartwig.hmftools.common.purple.repeat.RepeatContextFactory;
import com.hartwig.hmftools.common.region.GenomeRegion;
import com.hartwig.hmftools.common.region.GenomeRegionSelector;
import com.hartwig.hmftools.common.region.GenomeRegionSelectorFactory;
import com.hartwig.hmftools.common.variant.EnrichedSomaticVariant;
import com.hartwig.hmftools.common.variant.SomaticVariant;

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
    private final IndexedFastaSequenceFile reference;

    public EnrichedSomaticVariantFactory(double purity, double normFactor,
            @NotNull  final Multimap<String, GenomeRegion> highConfidenceRegions,
            @NotNull  final Multimap<String, PurpleCopyNumber> copyNumbers, @NotNull final IndexedFastaSequenceFile reference) {
        purityAdjuster = new PurityAdjuster(Gender.MALE, purity, normFactor);
        highConfidenceSelector = GenomeRegionSelectorFactory.create(highConfidenceRegions);
        copyNumberSelector = GenomeRegionSelectorFactory.create(copyNumbers);
        this.reference = reference;
    }

    public List<EnrichedSomaticVariant> enrich(final List<SomaticVariant> variants) {
        return variants.stream().map(this::enrich).collect(Collectors.toList());
    }

    private EnrichedSomaticVariant enrich(@NotNull final SomaticVariant variant) {
        final Builder builder = createBuilder(variant);

        highConfidenceSelector.select(variant).ifPresent(x -> inHighConfidenceRegion(builder));
        copyNumberSelector.select(variant).ifPresent(x -> addCopyNumber(builder, x, variant.alleleFrequency()));
        addTrinucleotideContext(builder, variant);
        addGenomeContext(builder, variant);

        return builder.build();
    }

    private Builder createBuilder(@NotNull final SomaticVariant variant) {
        return builder().from(variant)
                .trinucleotideContext("")
                .microhomology("")
                .refGenomeContext("")
                .repeatCount(0)
                .repeatSequence("")
                .totalReadCount(variant.totalReadCount())
                .alleleReadCount(variant.alleleReadCount())
                .highConfidenceRegion(false)
                .adjustedCopyNumber(0)
                .adjustedVAF(0);
    }

    private Builder addGenomeContext(@NotNull final Builder builder, @NotNull final SomaticVariant variant) {
        long positionBeforeEvent = variant.position();
        long start = Math.max(positionBeforeEvent - 100, 0);
        long end = positionBeforeEvent + 100;
        int position = (int) (positionBeforeEvent - start);
        final String sequence = reference.getSubsequenceAt(variant.chromosome(), start, end).getBaseString();
        builder.refGenomeContext(sequence);

        RepeatContextFactory.repeats(position, sequence, variant.ref(), variant.alt())
                .ifPresent(x -> builder.repeatSequence(x.sequence()).repeatCount(x.count()));

        // Deletions only atm
        if (variant.ref().length() > variant.alt().length()) {
            final String microhomology = Microhomology.microhomology(position, sequence, variant.ref(), variant.alt());
            return builder.microhomology(microhomology);
        }

        return builder;
    }

    private Builder addCopyNumber(@NotNull final Builder builder, @NotNull final PurpleCopyNumber copyNumber, double alleleFrequency) {
        double adjustedVAF = purityAdjuster.purityAdjustedVAF(Math.max(0.001, copyNumber.averageTumorCopyNumber()), alleleFrequency);
        return builder.adjustedCopyNumber(copyNumber.averageTumorCopyNumber()).adjustedVAF(adjustedVAF);
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
