package com.hartwig.hmftools.sage.quality;

import com.google.common.annotations.VisibleForTesting;
import com.hartwig.hmftools.sage.common.SimpleVariant;

import htsjdk.samtools.SAMRecord;

public class UltimaRealignedQualModel
{
    private final SimpleVariant mVariant;
    private final UltimaQualModel mBaseQualModel;
    private final int mVarReadIndexOffset;
    private final int mVarIndex;
    private final int mVariantRefIndex;

    private MsiJitterQualCache mMsiJitterQualCache;

    public UltimaRealignedQualModel(final SimpleVariant variant, final UltimaQualModel baseQualModel, int varReadIndexOffset, int varIndex, int variantRefIndex)
    {
        mVariant = variant;
        mBaseQualModel = baseQualModel;
        mVarReadIndexOffset = varReadIndexOffset;
        mVarIndex = varIndex;
        mVariantRefIndex = variantRefIndex;

        mMsiJitterQualCache = null;
    }

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

    @VisibleForTesting
    public UltimaRealignedQualModel(final SimpleVariant variant)
    {
        this(variant, -1);
    }

    public byte calculateQual(final SAMRecord record, final int varReadIndex)
    {
        return mBaseQualModel.calculateQual(record, varReadIndex + mVarReadIndexOffset);
    }

    public MsiJitterQualCache qualCache(final byte[] refBases, final byte[] readBases, final QualityCalculator qualityCalculator, final String sampleId)
    {
        if(mMsiJitterQualCache == null)
        {
            mMsiJitterQualCache = new MsiJitterQualCache(mVariant, mVarIndex, mVariantRefIndex, refBases, readBases, qualityCalculator, sampleId);
        }

        return mMsiJitterQualCache;
    }

    public UltimaQualModel baseQualModel() { return mBaseQualModel; }
    public int varReadIndexOffset() { return mVarReadIndexOffset; }
    public SimpleVariant variant() { return mVariant; }

    public UltimaModelType type()
    {
        // HOMOPOLYMER_ADJUSTMENT is a debugging placeholder
        return mBaseQualModel == null ? UltimaModelType.HOMOPOLYMER_ADJUSTMENT : mBaseQualModel.type();
    }
}