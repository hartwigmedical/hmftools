package com.hartwig.hmftools.bachelor.types;

import static com.hartwig.hmftools.bachelor.types.BachelorConfig.BACH_LOGGER;

import java.util.List;
import java.util.Optional;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.variant.Microhomology;
import com.hartwig.hmftools.common.variant.SomaticVariant;
import com.hartwig.hmftools.common.variant.VariantType;
import com.hartwig.hmftools.common.variant.repeat.RepeatContext;
import com.hartwig.hmftools.common.variant.repeat.RepeatContextFactory;

import org.apache.commons.math3.util.Pair;
import org.apache.logging.log4j.util.Strings;
import org.jetbrains.annotations.NotNull;

import htsjdk.samtools.reference.IndexedFastaSequenceFile;
import htsjdk.samtools.reference.ReferenceSequence;

public class RefGenomeEnrichedSomaticVariantFactory {

    private final IndexedFastaSequenceFile reference;

    RefGenomeEnrichedSomaticVariantFactory(@NotNull final IndexedFastaSequenceFile reference) {
        this.reference = reference;
    }

    @NotNull
    public List<EnrichedSomaticVariant> enrich(@NotNull List<? extends SomaticVariant> variants)
    {
        final List<EnrichedSomaticVariant> result = Lists.newArrayList();

        for (SomaticVariant variant : variants)
        {
            result.add(enrich(variant).build());
        }

        return result;
    }

    @NotNull
    protected ImmutableEnrichedSomaticVariant.Builder enrich(@NotNull final SomaticVariant variant)
    {
        final ImmutableEnrichedSomaticVariant.Builder builder = createBuilder(variant);
        addTrinucleotideContext(builder, variant, reference);
        addGenomeContext(builder, variant, reference);
        return builder;
    }

    @NotNull
    private static ImmutableEnrichedSomaticVariant.Builder createBuilder(@NotNull final SomaticVariant variant)
    {
        return ImmutableEnrichedSomaticVariant.builder().from(variant)
                .trinucleotideContext(Strings.EMPTY)
                .microhomology(Strings.EMPTY)
                .repeatCount(0)
                .repeatSequence(Strings.EMPTY)
                .highConfidenceRegion(false);
    }

    private static void addGenomeContext(
            @NotNull final ImmutableEnrichedSomaticVariant.Builder builder,
            @NotNull final SomaticVariant variant, @NotNull IndexedFastaSequenceFile reference)
    {
        final Pair<Integer, String> relativePositionAndRef = relativePositionAndRef(variant, reference);
        final Integer relativePosition = relativePositionAndRef.getFirst();
        final String sequence = relativePositionAndRef.getSecond();

        if (variant.type().equals(VariantType.INDEL))
        {
            if (variant.ref().length() > variant.alt().length())
            {
                builder.microhomology(Microhomology.microhomologyAtDelete(relativePosition, sequence, variant.ref()));
            }
            else if (variant.ref().length() == 1)
            {
                builder.microhomology(Microhomology.microhomologyAtInsert(relativePosition, sequence, variant.alt()));
            }
        }
        getRepeatContext(variant, relativePosition, sequence).ifPresent(x -> builder.repeatSequence(x.sequence()).repeatCount(x.count()));
    }

    @NotNull
    private static Pair<Integer, String> relativePositionAndRef(@NotNull final SomaticVariant variant,
            @NotNull final IndexedFastaSequenceFile reference)
    {
        final int chromosomeLength = reference.getSequenceDictionary().getSequence(variant.chromosome()).getSequenceLength();
        long positionBeforeEvent = variant.position();

        long start = Math.max(positionBeforeEvent - 100, 1);
        long end = Math.min(positionBeforeEvent + variant.ref().length() + 100 - 1, chromosomeLength - 1);
        int relativePosition = (int) (positionBeforeEvent - start);
        final String sequence;

        if (start < chromosomeLength && end < chromosomeLength)
        {
            sequence = reference.getSubsequenceAt(variant.chromosome(), start, end).getBaseString();
        }
        else
        {
            sequence = Strings.EMPTY;
            BACH_LOGGER.warn("Requested base sequence outside of chromosome region!");
        }

        return new Pair<>(relativePosition, sequence);
    }

    @NotNull
    private static Optional<RepeatContext> getRepeatContext(
            @NotNull final SomaticVariant variant, int relativePosition,
            @NotNull String sequence)
    {
        if (variant.type().equals(VariantType.INDEL))
        {
            return RepeatContextFactory.repeats(relativePosition + 1, sequence);
        }
        else if (variant.type().equals(VariantType.SNP) || variant.type().equals(VariantType.MNP))
        {
            Optional<RepeatContext> priorRepeat = RepeatContextFactory.repeats(relativePosition - 1, sequence);
            Optional<RepeatContext> postRepeat = RepeatContextFactory.repeats(relativePosition + variant.alt().length(), sequence);
            return max(priorRepeat, postRepeat);
        }
        else
        {
            return Optional.empty();
        }
    }

    @NotNull
    private static Optional<RepeatContext> max(
            @NotNull final Optional<RepeatContext> optionalPrior,
            @NotNull final Optional<RepeatContext> optionalPost)
    {
        if (!optionalPrior.isPresent())
            return optionalPost;

        if (!optionalPost.isPresent())
            return optionalPrior;

        final RepeatContext prior = optionalPrior.get();
        final RepeatContext post = optionalPost.get();

        if (post.sequence().length() > prior.sequence().length())
        {
            return optionalPost;
        }
        else if (post.sequence().length() == prior.sequence().length() && post.count() > prior.count())
        {
            return optionalPost;
        }

        return optionalPrior;
    }

    private static void addTrinucleotideContext(
            @NotNull final ImmutableEnrichedSomaticVariant.Builder builder,
            @NotNull final SomaticVariant variant, @NotNull IndexedFastaSequenceFile reference)
    {
        final int chromosomeLength = reference.getSequenceDictionary().getSequence(variant.chromosome()).getSequenceLength();

        if (variant.position() < chromosomeLength)
        {
            final ReferenceSequence sequence =
                    reference.getSubsequenceAt(variant.chromosome(), Math.max(1, variant.position() - 1), variant.position() + 1);
            builder.trinucleotideContext(sequence.getBaseString());
        }
        else
        {
            BACH_LOGGER.warn("Variant({}) position beyond ref genome", variant);
        }
    }
}
