package com.hartwig.hmftools.common.variant;

import java.util.Optional;
import java.util.function.Consumer;

import com.hartwig.hmftools.common.variant.enrich.VariantContextEnrichment;
import com.hartwig.hmftools.common.variant.repeat.RepeatContext;
import com.hartwig.hmftools.common.variant.repeat.RepeatContextFactory;

import org.apache.commons.math3.util.Pair;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.apache.logging.log4j.util.Strings;
import org.jetbrains.annotations.NotNull;

import htsjdk.samtools.reference.IndexedFastaSequenceFile;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.VariantContextBuilder;
import htsjdk.variant.vcf.VCFHeader;
import htsjdk.variant.vcf.VCFHeaderLineType;
import htsjdk.variant.vcf.VCFInfoHeaderLine;

public class RefContextEnrichment implements VariantContextEnrichment {

    private static final Logger LOGGER = LogManager.getLogger(RefContextEnrichment.class);

    public static final String TRINUCLEOTIDE_FLAG = "TRI";
    public static final String REPEAT_SEQUENCE_FLAG = "REPS";
    public static final String REPEAT_COUNT_FLAG = "REPC";
    public static final String MICROHOMOLOGY_FLAG = "MH";

    private static final String REPEAT_FLAG_DESCRIPTION = "Repeat sequence";
    private static final String TRI_FLAG_DESCRIPTION = "Trinucleotide context";
    private static final String MICROHOMOLOGY_FLAG_DESCRIPTION = "INDEL microhomology";
    private static final String REPEAT_COUNT_DESCRIPTION = "Repeat sequence count";

    private final IndexedFastaSequenceFile reference;
    private final Consumer<VariantContext> consumer;

    public RefContextEnrichment(@NotNull final IndexedFastaSequenceFile reference, final Consumer<VariantContext> consumer) {
        this.reference = reference;
        this.consumer = consumer;
    }

    @NotNull
    @Override
    public VCFHeader enrichHeader(@NotNull final VCFHeader template) {
        final VCFHeader outputVCFHeader = new VCFHeader(template.getMetaDataInInputOrder(), template.getSampleNamesInOrder());
        outputVCFHeader.addMetaDataLine(new VCFInfoHeaderLine(TRINUCLEOTIDE_FLAG, 1, VCFHeaderLineType.String, TRI_FLAG_DESCRIPTION));
        outputVCFHeader.addMetaDataLine(new VCFInfoHeaderLine(REPEAT_SEQUENCE_FLAG, 1, VCFHeaderLineType.String, REPEAT_FLAG_DESCRIPTION));
        outputVCFHeader.addMetaDataLine(new VCFInfoHeaderLine(REPEAT_COUNT_FLAG, 1, VCFHeaderLineType.Integer, REPEAT_COUNT_DESCRIPTION));
        outputVCFHeader.addMetaDataLine(new VCFInfoHeaderLine(MICROHOMOLOGY_FLAG,
                1,
                VCFHeaderLineType.String,
                MICROHOMOLOGY_FLAG_DESCRIPTION));

        return outputVCFHeader;
    }

    @Override
    public void accept(final VariantContext context) {
        final VariantContextBuilder builder = new VariantContextBuilder(context);
        final Pair<Integer, String> relativePositionAndRef = relativePositionAndRef(context);

        addTrinucleotideContext(builder, relativePositionAndRef);
        addMicrohomology(context, builder, relativePositionAndRef);
        addRepeatContext(context, builder, relativePositionAndRef);

        consumer.accept(builder.make());
    }

    @Override
    public void flush() {
        // None
    }

    private void addTrinucleotideContext(@NotNull final VariantContextBuilder builder,
            @NotNull final Pair<Integer, String> relativePositionAndRef) {
        final int relativePosition = relativePositionAndRef.getFirst();
        final String sequence = relativePositionAndRef.getSecond();
        if (!sequence.isEmpty()) {
            final String tri = sequence.substring(Math.max(0, relativePosition - 1), Math.min(sequence.length(), relativePosition + 2));
            builder.attribute(TRINUCLEOTIDE_FLAG, tri);
        }
    }

    private void addRepeatContext(@NotNull final VariantContext variant, @NotNull final VariantContextBuilder builder,
            final Pair<Integer, String> relativePositionAndRef) {

        final int relativePosition = relativePositionAndRef.getFirst();
        final String sequence = relativePositionAndRef.getSecond();

        Optional<RepeatContext> repeatContext = getRepeatContext(variant, relativePosition, sequence);
        if (repeatContext.isPresent()) {

            builder.attribute(REPEAT_SEQUENCE_FLAG, repeatContext.get().sequence());
            builder.attribute(REPEAT_COUNT_FLAG, repeatContext.get().count());
        }
    }

    private void addMicrohomology(@NotNull final VariantContext variant, @NotNull final VariantContextBuilder builder,
            final Pair<Integer, String> relativePositionAndRef) {
        final int relativePosition = relativePositionAndRef.getFirst();
        final String sequence = relativePositionAndRef.getSecond();
        if (variant.isIndel()) {
            final String ref = variant.getReference().getBaseString();
            final String alt = variant.getAlternateAllele(0).getBaseString();

            if (ref.length() > alt.length()) {
                builder.attribute(MICROHOMOLOGY_FLAG, Microhomology.microhomologyAtDelete(relativePosition, sequence, ref));
            } else if (ref.length() == 1) {
                builder.attribute(MICROHOMOLOGY_FLAG, Microhomology.microhomologyAtInsert(relativePosition, sequence, alt));
            }

        }
    }

    @NotNull
    private Pair<Integer, String> relativePositionAndRef(@NotNull final VariantContext variant) {
        final int refLength = variant.getReference().getBaseString().length();
        final int chromosomeLength = reference.getSequenceDictionary().getSequence(variant.getContig()).getSequenceLength();
        long positionBeforeEvent = variant.getStart();

        long start = Math.max(positionBeforeEvent - 100, 1);
        long end = Math.min(positionBeforeEvent + refLength + 100 - 1, chromosomeLength - 1);
        int relativePosition = (int) (positionBeforeEvent - start);
        final String sequence;
        if (start < chromosomeLength && end < chromosomeLength) {
            sequence = reference.getSubsequenceAt(variant.getContig(), start, end).getBaseString();
        } else {
            sequence = Strings.EMPTY;
            LOGGER.warn("Requested base sequence outside of chromosome region!");
        }
        return new Pair<>(relativePosition, sequence);
    }

    @NotNull
    private Optional<RepeatContext> getRepeatContext(@NotNull final VariantContext variant, int relativePosition,
            @NotNull final String sequence) {
        if (variant.isIndel()) {
            return RepeatContextFactory.repeats(relativePosition + 1, sequence);
        } else if (variant.isSNP() || variant.isMNP()) {
            int altLength = variant.getAlternateAllele(0).getBaseString().length();

            Optional<RepeatContext> priorRepeat = RepeatContextFactory.repeats(relativePosition - 1, sequence);
            Optional<RepeatContext> postRepeat = RepeatContextFactory.repeats(relativePosition + altLength, sequence);
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

}
