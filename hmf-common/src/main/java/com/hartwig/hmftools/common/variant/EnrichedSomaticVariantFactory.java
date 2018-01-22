package com.hartwig.hmftools.common.variant;

import static com.hartwig.hmftools.common.variant.ImmutableEnrichedSomaticVariant.Builder;
import static com.hartwig.hmftools.common.variant.ImmutableEnrichedSomaticVariant.builder;

import java.io.IOException;
import java.util.List;
import java.util.Optional;

import com.google.common.collect.Lists;
import com.google.common.collect.Multimap;
import com.hartwig.hmftools.common.purple.repeat.RepeatContext;
import com.hartwig.hmftools.common.purple.repeat.RepeatContextFactory;
import com.hartwig.hmftools.common.region.GenomeRegion;
import com.hartwig.hmftools.common.region.GenomeRegionSelector;
import com.hartwig.hmftools.common.region.GenomeRegionSelectorFactory;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.jetbrains.annotations.NotNull;

import htsjdk.samtools.reference.IndexedFastaSequenceFile;
import htsjdk.samtools.reference.ReferenceSequence;

public class EnrichedSomaticVariantFactory {

    private static final Logger LOGGER = LogManager.getLogger(EnrichedSomaticVariantFactory.class);

    @NotNull
    private final GenomeRegionSelector<GenomeRegion> highConfidenceSelector;
    @NotNull
    private final IndexedFastaSequenceFile reference;
    @NotNull
    private final ClonalityFactory clonalityFactory;

    private int unmatchedAnnotations;

    public EnrichedSomaticVariantFactory(@NotNull final Multimap<String, GenomeRegion> highConfidenceRegions,
            @NotNull final IndexedFastaSequenceFile reference, @NotNull final ClonalityFactory clonalityFactory) {
        this.highConfidenceSelector = GenomeRegionSelectorFactory.create(highConfidenceRegions);
        this.reference = reference;
        this.clonalityFactory = clonalityFactory;
    }

    @NotNull
    public List<EnrichedSomaticVariant> enrich(@NotNull List<PurityAdjustedSomaticVariant> variants) throws IOException {
        unmatchedAnnotations = 0;
        final List<EnrichedSomaticVariant> result = Lists.newArrayList();

        for (PurityAdjustedSomaticVariant variant : variants) {
            result.add(enrich(variant));
        }

        if (unmatchedAnnotations > 0) {
            LOGGER.warn("There were {} unmatched annotated genes.", unmatchedAnnotations);
        }

        return result;
    }

    @NotNull
    private EnrichedSomaticVariant enrich(@NotNull PurityAdjustedSomaticVariant variant) throws IOException {
        final Builder builder = createBuilder(variant);

        highConfidenceSelector.select(variant).ifPresent(x -> inHighConfidenceRegion(builder));
        addAnnotations(builder, variant);
        addTrinucleotideContext(builder, variant);
        addGenomeContext(builder, variant);
        builder.clonality(clonalityFactory.fromSample(variant));

        return builder.build();
    }

    private void addAnnotations(@NotNull final Builder builder, @NotNull final SomaticVariant variant) {
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
    }

    @NotNull
    private static Builder createBuilder(@NotNull final SomaticVariant variant) {
        return builder().from(variant)
                .trinucleotideContext("")
                .microhomology("")
                .gene("")
                .effect("")
                .repeatCount(0)
                .repeatSequence("")
                .highConfidenceRegion(false)
                .clonality(Clonality.UNKNOWN);
    }

    private void addGenomeContext(@NotNull final Builder builder, @NotNull final SomaticVariant variant) {
        long positionBeforeEvent = variant.position();
        long start = Math.max(positionBeforeEvent - 100, 1);
        long maxEnd = reference.getSequenceDictionary().getSequence(variant.chromosome()).getSequenceLength() - 1;
        long end = Math.min(positionBeforeEvent + 100, maxEnd);
        int relativePosition = (int) (positionBeforeEvent - start);
        final String sequence = reference.getSubsequenceAt(variant.chromosome(), start, end).getBaseString();
        if (variant.type().equals(VariantType.INDEL)) {
            builder.microhomology(Microhomology.microhomology(relativePosition, sequence, variant.ref(), variant.alt()));
        }
        getRepeatContext(variant, relativePosition, sequence).ifPresent(x -> builder.repeatSequence(x.sequence()).repeatCount(x.count()));
    }

    @NotNull
    public static Optional<RepeatContext> getRepeatContext(@NotNull final SomaticVariant variant, int relativePosition,
            @NotNull String sequence) {
        if (variant.type().equals(VariantType.INDEL)) {
            return RepeatContextFactory.repeats(relativePosition + 1, sequence);
        } else if (variant.type().equals(VariantType.SNP)) {
            Optional<RepeatContext> priorRepeat = RepeatContextFactory.repeats(relativePosition - 1, sequence);
            Optional<RepeatContext> postRepeat = RepeatContextFactory.repeats(relativePosition + 1, sequence);
            return max(priorRepeat, postRepeat);
        } else {
            return Optional.empty();
        }
    }

    @NotNull
    private static Optional<RepeatContext> max(@NotNull final Optional<RepeatContext> optionalPrior,
            @NotNull final Optional<RepeatContext> optionalPost) {
        if (!optionalPrior.isPresent()) {
            return optionalPost;
        }

        if (!optionalPost.isPresent()) {
            return optionalPrior;
        }

        final RepeatContext prior = optionalPrior.get();
        final RepeatContext post = optionalPost.get();

        if (post.sequence().length() > prior.sequence().length()) {
            return optionalPost;
        } else if (post.sequence().length() == prior.sequence().length() && post.count() > prior.count()) {
            return optionalPost;
        }

        return optionalPrior;
    }

    private void addTrinucleotideContext(@NotNull final Builder builder, @NotNull final SomaticVariant variant) {
        final ReferenceSequence sequence =
                reference.getSubsequenceAt(variant.chromosome(), Math.max(1, variant.position() - 1), variant.position() + 1);
        builder.trinucleotideContext(sequence.getBaseString());
    }

    @NotNull
    private static Builder inHighConfidenceRegion(@NotNull final Builder builder) {
        return builder.highConfidenceRegion(true);
    }
}
