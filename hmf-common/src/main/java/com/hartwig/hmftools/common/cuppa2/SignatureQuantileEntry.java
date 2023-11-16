package com.hartwig.hmftools.common.cuppa2;

public class SignatureQuantileEntry
{
    public final String FeatName;
    public final double FeatValue;
    public final String CancerType;
    public final double DataValue;
    public final int Rank;
    public final int RankGroup;

    SignatureQuantileEntry(
            final String featName,
            final double featValue,
            final String cancerType,
            final double dataValue,
            final int rank,
            final int rankGroup
    )
    {
        FeatName = featName;
        FeatValue = featValue;
        CancerType = cancerType;
        DataValue = dataValue;
        Rank = rank;
        RankGroup = rankGroup;
    }

    public static SignatureQuantileEntry fromCuppaPredictionEntry(CuppaPredictionEntry predictionEntry)
    {
        return new SignatureQuantileEntry(
                predictionEntry.FeatName,
                predictionEntry.FeatValue,
                predictionEntry.CancerType,
                predictionEntry.DataValue,
                predictionEntry.Rank,
                predictionEntry.RankGroup
        );
    }
}

