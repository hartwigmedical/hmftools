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
        if(taskIndex < mLastWrittenIndex)
            flushChromosome();

        if(taskIndex == mLastWrittenIndex + 1)
        {
            writeVariants(variants);
            mLastWrittenIndex = taskIndex;
            checkQueue();
        }
        else
        {
            int index = 0;
            while(index < mCompletedVariants.size())
            {
                if(taskIndex < mCompletedVariants.get(index).TaskIndex)
                    break;

                ++index;
            }

            mCompletedVariants.add(index, new CompleteVariants(taskIndex, variants));
        }
    }

    private void writeVariants(final List<SageVariant> variants)
    {
        variants.forEach(x -> mVcfFile.write(VariantContextFactory.create(x, mReferenceIds, mTumorIds)));
    }

    private void checkQueue()
    {
        int index = 0;

        while(index < mCompletedVariants.size())
        {
            CompleteVariants completeVariants = mCompletedVariants.get(index);

            if(completeVariants.TaskIndex > mLastWrittenIndex + 1)
                return;

            writeVariants(completeVariants.Variants);
            mLastWrittenIndex = completeVariants.TaskIndex;
            mCompletedVariants.remove(index);
        }
    }

    public void flushChromosome()
    {
        mCompletedVariants.forEach(x -> writeVariants(x.Variants));
        mCompletedVariants.clear();
        mLastWrittenIndex = -1;
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
