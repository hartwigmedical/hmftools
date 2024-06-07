package com.hartwig.hmftools.purple.germline;

import static com.hartwig.hmftools.common.variant.Hotspot.HOTSPOT_FLAG;
import static com.hartwig.hmftools.common.variant.CommonVcfTags.PASS;
import static com.hartwig.hmftools.common.variant.CommonVcfTags.REPORTED_FLAG;

import java.util.Set;

import com.google.common.collect.Sets;
import com.hartwig.hmftools.common.variant.VariantContextDecorator;
import com.hartwig.hmftools.common.variant.VariantType;
import com.hartwig.hmftools.common.variant.impact.VariantImpact;

import htsjdk.variant.variantcontext.Genotype;
import htsjdk.variant.variantcontext.VariantContext;

public class GermlineVariant
{
    private VariantContext mContext;

    private final String mChromosome;
    private VariantContextDecorator mDecorator;

    private final Set<String> mFilters; // modifiable

    public GermlineVariant(final VariantContext context)
    {
        mContext = context;
        mDecorator = new VariantContextDecorator(context);
        mChromosome = mContext.getContig();
        mFilters = Sets.newHashSet(context.getFilters());
    }

    public VariantContext context() { return mContext; }

    public void setContext(final VariantContext context)
    {
        mContext = context;
        mDecorator = new VariantContextDecorator(context);
    }

    public Set<String> filters() { return mFilters; }

    public String chromosome() { return mChromosome; }

    public final VariantContextDecorator decorator() { return mDecorator; }

    // convenience methods
    public VariantImpact variantImpact() { return mDecorator.variantImpact(); }

    public VariantType type() { return mDecorator.type(); }

    public boolean isPass() { return mFilters.isEmpty() || (mFilters.size() == 1 && mFilters.contains(PASS)); }

    public boolean isHotspot() { return mContext.hasAttribute(HOTSPOT_FLAG); }
    public String gene() { return mDecorator.variantImpact().GeneName; }
    public boolean biallelic() { return mDecorator.biallelic(); }
    public double copyNumber() { return mDecorator.variantCopyNumber(); }

    public boolean reported()
    {
        return mContext.hasAttribute(REPORTED_FLAG);
    }

    public Genotype getGenotype(final String sampleId)
    {
        return mContext.getGenotype(sampleId);
    }

}
