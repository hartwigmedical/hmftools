package com.hartwig.hmftools.common.basequal.jitter;

import com.hartwig.hmftools.common.qual.BqrReadType;

public class JitterModelParamsConsensus extends JitterModelParams
{
    public final BqrReadType ConsensusType;
    public JitterModelParamsConsensus(
            final String repeatUnit,
            final double optimalScaleRepeat4, final double optimalScaleRepeat5, final double optimalScaleRepeat6,
            final double scaleFitGradient, final double scaleFitIntercept, final double microsatelliteSkew,
            final BqrReadType consensusType)
    {
        super(repeatUnit, optimalScaleRepeat4, optimalScaleRepeat5, optimalScaleRepeat6, scaleFitGradient, scaleFitIntercept, microsatelliteSkew);
        ConsensusType = consensusType;
    }

    public JitterModelParamsConsensus(final JitterModelParams jitterModelParams, final BqrReadType consensusType)
    {
        super(jitterModelParams.RepeatUnit, jitterModelParams.OptimalScaleRepeat4, jitterModelParams.OptimalScaleRepeat5,
                jitterModelParams.OptimalScaleRepeat6, jitterModelParams.ScaleFitGradient, jitterModelParams.ScaleFitIntercept,
                jitterModelParams.MicrosatelliteSkew);
        ConsensusType = consensusType;
    }
}