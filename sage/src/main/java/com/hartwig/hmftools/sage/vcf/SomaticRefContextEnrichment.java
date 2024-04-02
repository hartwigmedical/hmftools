package com.hartwig.hmftools.sage.vcf;

import static com.hartwig.hmftools.common.variant.SageVcfTags.MICROHOMOLOGY_FLAG;
import static com.hartwig.hmftools.common.variant.SageVcfTags.REPEAT_COUNT_FLAG;
import static com.hartwig.hmftools.common.variant.SageVcfTags.REPEAT_SEQUENCE_FLAG;
import static com.hartwig.hmftools.common.variant.SageVcfTags.TRINUCLEOTIDE_FLAG;
import static com.hartwig.hmftools.common.variant.VariantUtils.relativePositionAndRef;

import java.util.Optional;
import java.util.function.Consumer;

import com.hartwig.hmftools.sage.common.Microhomology;
import com.hartwig.hmftools.common.variant.SageVcfTags;
import com.hartwig.hmftools.common.variant.enrich.VariantContextEnrichment;
import com.hartwig.hmftools.common.variant.repeat.RepeatContext;
import com.hartwig.hmftools.common.variant.repeat.RepeatContextFactory;

import org.apache.commons.math3.util.Pair;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.jetbrains.annotations.Nullable;

import htsjdk.samtools.reference.IndexedFastaSequenceFile;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.vcf.VCFHeader;

public class SomaticRefContextEnrichment implements VariantContextEnrichment
{
    private final IndexedFastaSequenceFile mRefGenome;

    @Nullable
    private final Consumer<VariantContext> mConsumer;

    public SomaticRefContextEnrichment(final IndexedFastaSequenceFile reference, final Consumer<VariantContext> consumer)
    {
        mRefGenome = reference;
        mConsumer = consumer;
    }

    @Override
    public VCFHeader enrichHeader(final VCFHeader template) { return SageVcfTags.addRefContextHeader(template); }

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
