package com.hartwig.hmftools.purple.somatic;

import static com.hartwig.hmftools.common.variant.Hotspot.HOTSPOT_FLAG;

import com.hartwig.hmftools.common.genome.position.GenomePosition;
import com.hartwig.hmftools.common.variant.AllelicDepth;
import com.hartwig.hmftools.common.variant.VariantContextDecorator;
import com.hartwig.hmftools.common.variant.VariantType;
import com.hartwig.hmftools.common.variant.impact.VariantImpact;

import htsjdk.variant.variantcontext.VariantContext;

public class SomaticVariant implements GenomePosition
{
    private VariantContext mContext;

    private final String mChromosome;
    private final int mPosition;
    private VariantContextDecorator mDecorator;
    private final AllelicDepth mTumorAllelicDepth;
    private final AllelicDepth mReferenceAllelicDepth;

    public SomaticVariant(final VariantContext context, final String sampleId, final String referenceId)
    {
        mContext = context;
        mDecorator = new VariantContextDecorator(mContext);

        mChromosome = mContext.getContig();
        mPosition = mContext.getStart();
        mTumorAllelicDepth = sampleId != null ? mDecorator.allelicDepth(sampleId) :  null;
        mReferenceAllelicDepth = referenceId != null ? mDecorator.allelicDepth(referenceId) :  null;
    }

    public VariantContext context() { return mContext; }

    public void setContext(final VariantContext context)
    {
        mContext = context;
        mDecorator = new VariantContextDecorator(mContext);
    }

    @Override
    public String chromosome() { return mChromosome; }

    @Override
    public int position() { return mPosition; }

    // convenience methods
    public VariantContextDecorator decorator() { return mDecorator; }
    public VariantImpact variantImpact() { return mDecorator.variantImpact(); }

    public VariantType type() { return mDecorator.type(); }

    public boolean isPass() { return mDecorator.isPass(); }
    public boolean isFiltered() { return !isPass(); }

    public double copyNumber() { return mDecorator.variantCopyNumber(); }

    public boolean isHotspot() { return mContext.hasAttribute(HOTSPOT_FLAG); }
    public boolean biallelic() { return mDecorator.biallelic(); }
    public String gene() { return mDecorator.variantImpact().GeneName; }

    public boolean hasTumorAlleleDepth() { return mTumorAllelicDepth != null; }
    public AllelicDepth tumorAlleleDepth() { return mTumorAllelicDepth; }
    public AllelicDepth referenceAlleleDepth() { return mReferenceAllelicDepth; }
    public double alleleFrequency() { return mTumorAllelicDepth != null ? mTumorAllelicDepth.alleleFrequency() : 0; }
    public int totalReadCount() { return mTumorAllelicDepth != null ? mTumorAllelicDepth.TotalReadCount : 0; }
    public int alleleReadCount() { return mTumorAllelicDepth != null ? mTumorAllelicDepth.AlleleReadCount : 0; }
    public int referenceAlleleReadCount() { return mReferenceAllelicDepth != null ? mReferenceAllelicDepth.AlleleReadCount : 0; }

    public String toString()
    {
        return String.format("%s %s:%d %s>%s filter(%s) tier(%s) codingEffect(%s)",
                type(), chromosome(), position(), mDecorator.ref(), mDecorator.alt(), mDecorator.filter(), mDecorator.tier(),
                mDecorator.variantImpact().CanonicalCodingEffect);
    }
}
