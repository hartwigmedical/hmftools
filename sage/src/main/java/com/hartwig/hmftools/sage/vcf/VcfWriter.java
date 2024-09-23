package com.hartwig.hmftools.sage.vcf;

import java.util.List;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.sage.SageConfig;
import com.hartwig.hmftools.sage.common.SageVariant;

import htsjdk.samtools.reference.IndexedFastaSequenceFile;

public class VcfWriter
{
    private final List<String> mTumorIds;
    private final List<String> mReferenceIds;
    private final VariantVCF mVcfFile;

    // state to write variants in order
    private int mLastWrittenIndex;
    private final List<CompleteVariants> mCompletedVariants;

    public VcfWriter(
            final SageConfig config, final List<String> tumorIds, final List<String> referenceIds, final IndexedFastaSequenceFile refGenome)
    {
        mTumorIds = tumorIds;
        mReferenceIds = referenceIds;
        mVcfFile = new VariantVCF(refGenome, config, tumorIds, referenceIds);
        mCompletedVariants = Lists.newArrayList();
        mLastWrittenIndex = -1;
    }

    public void writeVariants(int taskIndex, final List<SageVariant> variants)
    {
        // queue the newly completed variants then check which completed regions can be written
        int index = 0;
        while(index < mCompletedVariants.size())
        {
            if(taskIndex < mCompletedVariants.get(index).TaskIndex)
                break;

            ++index;
        }

        mCompletedVariants.add(index, new CompleteVariants(taskIndex, variants));

        checkQueue();
    }

    private void checkQueue()
    {
        // write up until the set of second last successive completed region
        // eg if regions 0, 1, 3 are queued then only 0 will be written, and then 1 will be written when 2 arrives
        // this is so variants in a next region with earlier spanning positions can be included in the one being written
        int index = 0;
        while(index < mCompletedVariants.size() - 1)
        {
            CompleteVariants completeVariants = mCompletedVariants.get(index);

            if(completeVariants.TaskIndex > mLastWrittenIndex + 1) // break if there is a gap
                return;

            // break if the next set of variants hasn't completed yet
            CompleteVariants nextCompleteVariants = mCompletedVariants.get(index + 1);

            if(nextCompleteVariants.TaskIndex > completeVariants.TaskIndex + 1)
                break;

            checkLaterVariants(completeVariants, nextCompleteVariants);

            writeVariants(completeVariants.Variants);
            mLastWrittenIndex = completeVariants.TaskIndex;
            mCompletedVariants.remove(index);
        }
    }

    private void checkLaterVariants(final CompleteVariants completeVariants, final CompleteVariants nextCompleteVariants)
    {
        // some variants at the start of a region could have their position adjusted to fall into any earlier region
        while(!nextCompleteVariants.Variants.isEmpty())
        {
            SageVariant variant = nextCompleteVariants.Variants.get(0);

            int lastIndex = completeVariants.Variants.size() - 1;
            if(variant.position() >= completeVariants.Variants.get(lastIndex).position())
                break;

            // merge in this later variant into the earlier completed set of variants into its correct position
            nextCompleteVariants.Variants.remove(0);

            int index = lastIndex;
            while(index >= 0)
            {
                SageVariant earlierVariant = completeVariants.Variants.get(index);

                if(earlierVariant.position() <= variant.position())
                    break;

                --index;
            }

            completeVariants.Variants.add(index + 1, variant);
        }
    }

    public void flushChromosome()
    {
        mCompletedVariants.forEach(x -> writeVariants(x.Variants));
        mCompletedVariants.clear();
        mLastWrittenIndex = -1;
    }

    private void writeVariants(final List<SageVariant> variants)
    {
        variants.forEach(x -> mVcfFile.write(VariantContextFactory.create(x, mReferenceIds, mTumorIds)));
    }

    public void close()
    {
        mVcfFile.close();
    }

    private class CompleteVariants
    {
        public final int TaskIndex;
        public final List<SageVariant> Variants;

        public CompleteVariants(int taskIndex, final List<SageVariant> variants)
        {
            Variants = variants;
            TaskIndex = taskIndex;
        }

        public String toString() { return String.format("%d: %d variants", TaskIndex, Variants.size()); }
    }
}
