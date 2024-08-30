// TODO: REVIEW
package com.hartwig.hmftools.sage.quality;

import static com.hartwig.hmftools.sage.quality.QualityCalculator.INVALID_BASE_QUAL;

import com.hartwig.hmftools.common.qual.BaseQualAdjustment;
import com.hartwig.hmftools.sage.common.SimpleVariant;
import com.hartwig.hmftools.sage.common.VariantReadContext;

public class MsiJitterQualCache
{
    private final double mMsiIndelErrorQual;
    private final boolean mIsMsiSampleAndVariant;

    public MsiJitterQualCache(final VariantReadContext readContext, final QualityCalculator qualityCalculator, final String sampleId)
    {
        this(readContext.variant(), readContext.VarIndex, readContext.variantRefIndex(), readContext.RefBases, readContext.ReadBases, qualityCalculator, sampleId);
    }

    public MsiJitterQualCache(final SimpleVariant variant, int varIndex, int variantRefIndex, final byte[] refBases, final byte[] readBases,
            final QualityCalculator qualityCalculator, final String sampleId)
    {
        double errorRate = qualityCalculator.msiJitterCalcs().calcErrorRate(
                variant, varIndex, variantRefIndex, refBases, readBases, sampleId);

        mMsiIndelErrorQual = errorRate > 0 ? BaseQualAdjustment.probabilityToPhredQual(errorRate) : INVALID_BASE_QUAL;
        mIsMsiSampleAndVariant = usesMsiIndelErrorQual() && qualityCalculator.msiJitterCalcs().getProbableMsiStatus(sampleId);
    }

    public double msiIndelErrorQual() { return mMsiIndelErrorQual; }
    public boolean usesMsiIndelErrorQual() { return mMsiIndelErrorQual != INVALID_BASE_QUAL; }
    public boolean isMsiSampleAndVariant() { return mIsMsiSampleAndVariant; }
}
