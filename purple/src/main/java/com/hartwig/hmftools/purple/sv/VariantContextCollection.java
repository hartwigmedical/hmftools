package com.hartwig.hmftools.purple.sv;

import java.util.Collections;
import java.util.Iterator;
import java.util.List;
import java.util.TreeSet;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.sv.StructuralVariant;
import com.hartwig.hmftools.common.sv.StructuralVariantFactory;

import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.VariantContextComparator;
import htsjdk.variant.vcf.VCFHeader;

public class VariantContextCollection
{
    private final TreeSet<VariantContext> mVariantContexts;
    private final List<StructuralVariant> mVariants;

    private boolean mModified;

    public VariantContextCollection(final VCFHeader header)
    {
        if(header != null)
            mVariantContexts = new TreeSet<>(new VCComparator(header.getSequenceDictionary()));
        else
            mVariantContexts = null;

        mVariants = Lists.newArrayList();
    }

    public void clear()
    {
        if(mVariantContexts != null)
            mVariantContexts.clear();

        mVariants.clear();
    }

    public void addVariant(final VariantContext variantContext)
    {
        if(mVariantContexts == null)
            return;

        mModified = true;
        mVariantContexts.add(variantContext);
    }

    public List<StructuralVariant> variants()
    {
        if(mVariantContexts == null)
            return Collections.emptyList();

        if(mModified)
        {
            // converts variant contexts into structural variants
            mModified = false;
            StructuralVariantFactory svFactory = StructuralVariantFactory.build(new SegmentationVariantsFilter());
            mVariantContexts.forEach(svFactory::addVariantContext);

            mVariants.clear();
            mVariants.addAll(svFactory.results());
        }

        return mVariants;
    }

    public Iterator<VariantContext> iterator()
    {
        return mVariantContexts != null ? mVariantContexts.iterator() : Collections.emptyIterator();
    }

    private static class VCComparator extends VariantContextComparator
    {
        VCComparator(final SAMSequenceDictionary dictionary)
        {
            super(dictionary);
        }

        @Override
        public int compare(final VariantContext firstVariantContext, final VariantContext secondVariantContext)
        {
            int positionResult = super.compare(firstVariantContext, secondVariantContext);

            return positionResult == 0 ? firstVariantContext.getID().compareTo(secondVariantContext.getID()) : positionResult;
        }
    }
}
