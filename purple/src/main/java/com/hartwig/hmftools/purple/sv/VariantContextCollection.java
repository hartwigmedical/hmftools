package com.hartwig.hmftools.purple.sv;

import static com.hartwig.hmftools.common.sv.SvFactoryInterface.buildSvFactory;

import java.util.Collections;
import java.util.Iterator;
import java.util.List;
import java.util.TreeSet;
import java.util.function.Predicate;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.sv.StructuralVariant;
import com.hartwig.hmftools.common.sv.SvFactoryInterface;

import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.VariantContextComparator;
import htsjdk.variant.vcf.VCFHeader;

public class VariantContextCollection
{
    private final TreeSet<VariantContext> mVariantContexts;
    private final List<StructuralVariant> mVariants;
    private final boolean mUseGridssSVs;

    private boolean mModified;

    public VariantContextCollection(final VCFHeader header, boolean useGridssSVs)
    {
        mUseGridssSVs = useGridssSVs;

        if(header != null)
            mVariantContexts = new TreeSet<>(new VCComparator(header.getSequenceDictionary()));
        else
            mVariantContexts = null;

        mVariants = Lists.newArrayList();
    }

    public void add(final VariantContext variantContext)
    {
        if(mVariantContexts == null)
            return;

        mModified = true;
        if(variantContext.contains(variantContext))
        {
            mVariantContexts.remove(variantContext);
            mVariantContexts.add(variantContext);
        }
    }

    public int remove(final Predicate<VariantContext> removePredicate)
    {
        if(mVariantContexts == null)
            return 0;

        int removed = 0;
        final Iterator<VariantContext> iterator = mVariantContexts.iterator();
        while(iterator.hasNext())
        {
            final VariantContext variantContext = iterator.next();
            if(removePredicate.test(variantContext))
            {
                iterator.remove();
                removed++;
                mModified = true;
            }
        }
        return removed;
    }

    public List<StructuralVariant> variants()
    {
        if(mVariantContexts == null)
            return Collections.emptyList();

        if(mModified)
        {
            // converts variant contexts into structural variants
            mModified = false;
            SvFactoryInterface svFactory = buildSvFactory(mUseGridssSVs, new SegmentationVariantsFilter());
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
