package com.hartwig.hmftools.svanalysis.types;

import com.hartwig.hmftools.common.variant.structural.annotation.GeneFusion;

public class SvFusion
{
    private RnaFusionData mRnaFusionData;
    private GeneFusion mFusion;
    private SvBreakend mUpBreakend;
    private SvBreakend mDownBreakend;

    public SvFusion(final GeneFusion fusion, final SvBreakend upBreakend, final SvBreakend downBreakend)
    {
        mFusion = fusion;
        mUpBreakend = upBreakend;
        mDownBreakend = downBreakend;

        mRnaFusionData = null;
    }

    public void setRnaFusionData(final RnaFusionData rnaFusionData)
    {
        mRnaFusionData = rnaFusionData;

    }

    public final RnaFusionData getRnaFusionData() { return mRnaFusionData; }
}
