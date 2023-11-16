package com.hartwig.hmftools.common.cuppa2;

public class FeatureContributionEntry
{
    public final String ClfName;
    public final String FeatName;
    public final double FeatValue;
    public final String CancerType;
    public final double DataValue;
    public final int Rank;
    public final int RankGroup;

    FeatureContributionEntry(
            final String clfName,
            final String featName,
            final double featValue,
            final String cancerType,
            final double dataValue,
            final int rank,
            final int rankGroup
    )
    {
        ClfName = clfName;
        FeatName = featName;
        FeatValue = featValue;
        CancerType = cancerType;
        DataValue = dataValue;
        Rank = rank;
        RankGroup = rankGroup;
    }

    public static FeatureContributionEntry fromCuppaPredictionEntry(CuppaPredictionEntry predictionEntry)
    {
        return new FeatureContributionEntry(
                predictionEntry.ClfName.toString(),
                predictionEntry.FeatName,
                predictionEntry.FeatValue,
                predictionEntry.CancerType,
                predictionEntry.DataValue,
                predictionEntry.Rank,
                predictionEntry.RankGroup
        );
    }
}

