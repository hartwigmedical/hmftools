package com.hartwig.hmftools.sage.seqtech;

import static java.lang.String.format;

import com.hartwig.hmftools.common.variant.SimpleVariant;
import com.hartwig.hmftools.sage.quality.MsiJitterQualCache;
import com.hartwig.hmftools.sage.quality.QualityCalculator;

import htsjdk.samtools.SAMRecord;

public class UltimaRealignedQualModel
{
    private final SimpleVariant mVariant;
    private final UltimaQualModel mQualModel;
    private final int mVarReadIndexOffset;
    private final int mVarIndex;
    private final int mVariantRefIndex;

    private MsiJitterQualCache mMsiJitterQualCache;

    public UltimaRealignedQualModel(
            final SimpleVariant variant, final UltimaQualModel baseQualModel, int varReadIndexOffset, int varIndex, int variantRefIndex)
    {
        mVariant = variant;
        mQualModel = baseQualModel;
        mVarReadIndexOffset = varReadIndexOffset;
        mVarIndex = varIndex;
        mVariantRefIndex = variantRefIndex;

        mMsiJitterQualCache = null;
    }

    public UltimaQualModel qualModel() { return mQualModel; }

    /*
    private UltimaRealignedQualModel(final SimpleVariant variant, int varReadIndexOffset, int varIndex, int variantRefIndex)
    {
        mVariant = variant;
        mBaseQualModel = null;
        mVarReadIndexOffset = varReadIndexOffset;
        mVarIndex = varIndex;
        mVariantRefIndex = variantRefIndex;

        mMsiJitterQualCache = null;
    }

    @VisibleForTesting
    public UltimaRealignedQualModel(final SimpleVariant variant, int varReadIndexOffset)
    {
        this(variant, varReadIndexOffset, -1, -1);
    }
    */

    public byte calculateQual(final SAMRecord record, final int varReadIndex)
    {
        return mQualModel.calculateQual(record, varReadIndex + mVarReadIndexOffset);
    }

    public MsiJitterQualCache qualCache(
            final byte[] refBases, final byte[] readBases, final QualityCalculator qualityCalculator, final String sampleId)
    {
        if(mMsiJitterQualCache == null)
        {
            mMsiJitterQualCache = new MsiJitterQualCache(
                    mVariant, mVarIndex, mVariantRefIndex, refBases, readBases, qualityCalculator, sampleId);
        }

        return mMsiJitterQualCache;
    }

    public UltimaQualModel baseQualModel() { return mQualModel; }
    public int varReadIndexOffset() { return mVarReadIndexOffset; }
    public SimpleVariant variant() { return mVariant; }

    public UltimaModelType type() { return mQualModel.type(); } // mBaseQualModel == null ? UltimaModelType.NONE :

    public String toString()
    {
        return format("%s model(%s) indices(var=%d refIndex=%d readIndexOffset=%d)",
                mVariant, type(), mVarIndex, mVariantRefIndex, mVarReadIndexOffset);
    }
}