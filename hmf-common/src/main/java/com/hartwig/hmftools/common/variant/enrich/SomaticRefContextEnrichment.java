package com.hartwig.hmftools.common.variant.enrich;

import java.util.Optional;
import java.util.function.Consumer;

import com.hartwig.hmftools.common.variant.Microhomology;
import com.hartwig.hmftools.common.variant.repeat.RepeatContext;
import com.hartwig.hmftools.common.variant.repeat.RepeatContextFactory;

import org.apache.commons.math3.util.Pair;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.apache.logging.log4j.util.Strings;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

import htsjdk.samtools.SAMSequenceRecord;
import htsjdk.samtools.reference.IndexedFastaSequenceFile;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.vcf.VCFHeader;
import htsjdk.variant.vcf.VCFHeaderLineType;
import htsjdk.variant.vcf.VCFInfoHeaderLine;

public class SomaticRefContextEnrichment implements VariantContextEnrichment
{
    private static final Logger LOGGER = LogManager.getLogger(SomaticRefContextEnrichment.class);

    public static final String MICROHOMOLOGY_FLAG = "MH";
    public static final String TRINUCLEOTIDE_FLAG = "TNC";
    public static final String REPEAT_COUNT_FLAG = "REP_C";
    public static final String REPEAT_SEQUENCE_FLAG = "REP_S";

    private static final String REPEAT_FLAG_DESCRIPTION = "Repeat sequence";
    private static final String MICROHOMOLOGY_FLAG_DESCRIPTION = "Microhomology";
    private static final String REPEAT_COUNT_DESCRIPTION = "Repeat sequence count";
    private static final String TRINUCLEOTIDE_FLAG_DESCRIPTION = "Tri-nucleotide context";

    private final IndexedFastaSequenceFile mRefGenome;

    @Nullable
    private final Consumer<VariantContext> mConsumer;

    public SomaticRefContextEnrichment(final IndexedFastaSequenceFile reference, final Consumer<VariantContext> consumer)
    {
        mRefGenome = reference;
        mConsumer = consumer;
    }

    @Override
    public VCFHeader enrichHeader(final VCFHeader template) { return addHeader(template); }

    public static VCFHeader addHeader(final VCFHeader template)
    {
        template.addMetaDataLine(new VCFInfoHeaderLine(TRINUCLEOTIDE_FLAG, 1, VCFHeaderLineType.String, TRINUCLEOTIDE_FLAG_DESCRIPTION));
        template.addMetaDataLine(new VCFInfoHeaderLine(REPEAT_SEQUENCE_FLAG, 1, VCFHeaderLineType.String, REPEAT_FLAG_DESCRIPTION));
        template.addMetaDataLine(new VCFInfoHeaderLine(REPEAT_COUNT_FLAG, 1, VCFHeaderLineType.Integer, REPEAT_COUNT_DESCRIPTION));
        template.addMetaDataLine(new VCFInfoHeaderLine(MICROHOMOLOGY_FLAG, 1, VCFHeaderLineType.String, MICROHOMOLOGY_FLAG_DESCRIPTION));

        return template;
    }

    @Override
    public void accept(final VariantContext context)
    {
        processVariant(context);

        if(mConsumer != null)
            mConsumer.accept(context);
    }

    public void processVariant(final VariantContext context)
    {
        final Pair<Integer, String> relativePositionAndRef = relativePositionAndRef(mRefGenome, context);
        if(relativePositionAndRef.getFirst() > -1)
        {
            addTrinucleotideContext(context, relativePositionAndRef);
            addMicrohomology(context, relativePositionAndRef);
            addRepeatContext(context, relativePositionAndRef);
        }
    }

    @Override
    public void flush() { }

    private void addTrinucleotideContext(final VariantContext variant, final Pair<Integer, String> relativePositionAndRef)
    {
        final int relativePosition = relativePositionAndRef.getFirst();
        final String sequence = relativePositionAndRef.getSecond();
        if(!sequence.isEmpty())
        {
            final String tri = sequence.substring(Math.max(0, relativePosition - 1), Math.min(sequence.length(), relativePosition + 2));
            variant.getCommonInfo().putAttribute(TRINUCLEOTIDE_FLAG, tri, true);
        }
    }

    private void addRepeatContext(final VariantContext variant, final Pair<Integer, String> relativePositionAndRef)
    {
        final int relativePosition = relativePositionAndRef.getFirst();
        final String sequence = relativePositionAndRef.getSecond();

        Optional<RepeatContext> repeatContext = getRepeatContext(variant, relativePosition, sequence);
        if(repeatContext.isPresent())
        {
            variant.getCommonInfo().putAttribute(REPEAT_SEQUENCE_FLAG, repeatContext.get().sequence(), true);
            variant.getCommonInfo().putAttribute(REPEAT_COUNT_FLAG, repeatContext.get().count(), true);
        }
    }

    private void addMicrohomology(final VariantContext variant, final Pair<Integer, String> relativePositionAndRef)
    {
        final int relativePosition = relativePositionAndRef.getFirst();
        final String sequence = relativePositionAndRef.getSecond();
        if(variant.isIndel())
        {
            final String ref = variant.getReference().getBaseString();
            final String alt = variant.getAlternateAllele(0).getBaseString();

            final String microhomology = ref.length() > alt.length()
                    ? Microhomology.microhomologyAtDelete(relativePosition, sequence, ref)
                    : Microhomology.microhomologyAtInsert(relativePosition, sequence, alt);

            if(!microhomology.isEmpty())
            {
                variant.getCommonInfo().putAttribute(MICROHOMOLOGY_FLAG, microhomology, true);
            }
        }
    }

    public static Pair<Integer, String> relativePositionAndRef(final IndexedFastaSequenceFile reference, final VariantContext variant)
    {
        final int refLength = variant.getReference().getBaseString().length();
        final SAMSequenceRecord samSequenceRecord = reference.getSequenceDictionary().getSequence(variant.getContig());
        if(samSequenceRecord == null)
        {
            return new Pair<>(-1, Strings.EMPTY);
        }

        final int chromosomeLength = samSequenceRecord.getSequenceLength();
        int positionBeforeEvent = variant.getStart();

        int start = Math.max(positionBeforeEvent - 100, 1);
        int end = Math.min(positionBeforeEvent + refLength + 100 - 1, chromosomeLength - 1);
        int relativePosition = positionBeforeEvent - start;
        final String sequence;
        if(start < chromosomeLength && end < chromosomeLength)
        {
            sequence = reference.getSubsequenceAt(variant.getContig(), start, end).getBaseString();
        }
        else
        {
            sequence = Strings.EMPTY;
            LOGGER.warn("Requested base sequence outside of chromosome region!");
        }
        return new Pair<>(relativePosition, sequence);
    }

    private Optional<RepeatContext> getRepeatContext(final VariantContext variant, int relativePosition, final String sequence)
    {
        if(variant.isIndel())
        {
            return RepeatContextFactory.repeats(relativePosition + 1, sequence);
        }
        else if(variant.isSNP() || variant.isMNP())
        {
            int altLength = variant.getAlternateAllele(0).getBaseString().length();

            Optional<RepeatContext> priorRepeat = RepeatContextFactory.repeats(relativePosition - 1, sequence);
            Optional<RepeatContext> postRepeat = RepeatContextFactory.repeats(relativePosition + altLength, sequence);
            return max(priorRepeat, postRepeat);
        }
        else
        {
            return Optional.empty();
        }
    }

    private static Optional<RepeatContext> max(final Optional<RepeatContext> optionalPrior, final Optional<RepeatContext> optionalPost)
    {
        if(!optionalPrior.isPresent())
        {
            return optionalPost;
        }

        if(!optionalPost.isPresent())
        {
            return optionalPrior;
        }

        final RepeatContext prior = optionalPrior.get();
        final RepeatContext post = optionalPost.get();

        if(post.sequence().length() > prior.sequence().length())
        {
            return optionalPost;
        }
        else if(post.sequence().length() == prior.sequence().length() && post.count() > prior.count())
        {
            return optionalPost;
        }

        return optionalPrior;
    }
}
