package com.hartwig.hmftools.common.cuppa2;

public class ProbabilityEntry
{
    public final String ClfName;
    public final String CancerType;
    public final double DataValue;
    public final int Rank;
    public final int RankGroup;

    ProbabilityEntry(
            final String clfName,
            final String cancerType,
            final double dataValue,
            final int rank,
            final int rankGroup
    )
    {
        ClfName = clfName;
        CancerType = cancerType;
        DataValue = dataValue;
        Rank = rank;
        RankGroup = rankGroup;
    }

    public static ProbabilityEntry fromCuppaPredictionEntry(CuppaPredictionEntry predictionEntry)
    {
        return new ProbabilityEntry(
                predictionEntry.ClfName.toString(),
                predictionEntry.CancerType,
                predictionEntry.DataValue,
                predictionEntry.Rank,
                predictionEntry.RankGroup
        );
    }
}

