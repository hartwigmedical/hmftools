package com.hartwig.hmftools.bachelor;

import static com.hartwig.hmftools.bachelor.types.BachelorConfig.BACH_LOGGER;

import java.util.Optional;

import com.hartwig.hmftools.common.variant.Microhomology;
import com.hartwig.hmftools.common.variant.VariantType;
import com.hartwig.hmftools.common.variant.repeat.RepeatContext;
import com.hartwig.hmftools.common.variant.repeat.RepeatContextFactory;

import org.jetbrains.annotations.NotNull;

import htsjdk.samtools.reference.IndexedFastaSequenceFile;
import htsjdk.samtools.reference.ReferenceSequence;

public class RefGenomeEnrichment
{
    private final IndexedFastaSequenceFile mRefGenome;

    public RefGenomeEnrichment(final IndexedFastaSequenceFile refGenome)
    {
        mRefGenome = refGenome;
    }

    public boolean isValid() { return mRefGenome != null; }

    public String getMicrohomology(final VariantType variantType, final String ref, final String alt, long position, final String sequence)
    {
        if(variantType.equals(VariantType.INDEL))
        {
            int relativePosition = getRelativePosition(position);

            if(ref.length() > alt.length())
            {
                return Microhomology.microhomologyAtDelete(relativePosition, sequence, ref);
            }
            else if(ref.length() == 1)
            {
                return Microhomology.microhomologyAtInsert(relativePosition, sequence, alt);
            }
        }

        return "";
    }

    private int getRelativePosition(final long position)
    {
        long start = Math.max(position - 100, 1);
        return (int)(position - start);
    }

    public final String getSequence(final String chromosome, final long position, final String ref)
    {
        int chromosomeLength = mRefGenome.getSequenceDictionary().getSequence(chromosome).getSequenceLength();

        long start = Math.max(position - 100, 1);
        long end = Math.min(position + ref.length() + 100 - 1, chromosomeLength - 1);

        if (start < chromosomeLength && end < chromosomeLength)
            return mRefGenome.getSubsequenceAt(chromosome, start, end).getBaseString();

        BACH_LOGGER.warn("chr({}) position({}) outside genome", chromosome, position);
        return "";
    }

    public Optional<RepeatContext> getRepeatContext(final VariantType variantType, final String sequence, long position, final String alt)
    {
        final int relativePosition = getRelativePosition(position);

        if (variantType.equals(VariantType.INDEL))
        {
            return RepeatContextFactory.repeats(relativePosition + 1, sequence);
        }
        else if (variantType.equals(VariantType.SNP) || variantType.equals(VariantType.MNP))
        {
            Optional<RepeatContext> priorRepeat = RepeatContextFactory.repeats(relativePosition - 1, sequence);
            Optional<RepeatContext> postRepeat = RepeatContextFactory.repeats(relativePosition + alt.length(), sequence);
            return max(priorRepeat, postRepeat);
        }
        else
        {
            return Optional.empty();
        }
    }

    private static Optional<RepeatContext> max(
            @NotNull final Optional<RepeatContext> optionalPrior, @NotNull final Optional<RepeatContext> optionalPost)
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

    public String getTrinucleotideContext(final String chromosome, long position)
    {
        final int chromosomeLength = mRefGenome.getSequenceDictionary().getSequence(chromosome).getSequenceLength();

        if (position < chromosomeLength)
        {
            final ReferenceSequence sequence =
                    mRefGenome.getSubsequenceAt(chromosome, Math.max(1, position - 1), position + 1);

            return sequence.getBaseString();
        }
        else
        {
            BACH_LOGGER.warn("Variant chr({}) position({}) beyond ref genome", chromosome, position);
        }

        return "";
    }

}
