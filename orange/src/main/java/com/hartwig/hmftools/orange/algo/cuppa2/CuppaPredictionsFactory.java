package com.hartwig.hmftools.orange.algo.cuppa2;

import java.util.ArrayList;
import java.util.List;

import com.hartwig.hmftools.common.cuppa2.Categories;
import com.hartwig.hmftools.common.cuppa2.CuppaPredictionEntry;
import com.hartwig.hmftools.common.cuppa2.CuppaPredictions;
import com.hartwig.hmftools.datamodel.cuppa2.ImmutableProbabilityEntry;
import com.hartwig.hmftools.datamodel.cuppa2.ImmutableFeatureContributionEntry;
import com.hartwig.hmftools.datamodel.cuppa2.ImmutableSignatureQuantileEntry;
import com.hartwig.hmftools.datamodel.cuppa2.ImmutableCuppaPredictions;

import org.jetbrains.annotations.NotNull;

public class CuppaPredictionsFactory
{
    @NotNull
    public static ImmutableCuppaPredictions create(@NotNull CuppaPredictions cuppaPredictions)
    {
        return ImmutableCuppaPredictions.builder()
                .topPrediction(getTopPrediction(cuppaPredictions))
                .probs(getProbabilities(cuppaPredictions))
                .featContribs(getFeatureContributions(cuppaPredictions))
                .sigQuantiles(getSignatureQuantiles(cuppaPredictions))
                .build();
    }

    @NotNull
    public static ImmutableProbabilityEntry getTopPrediction(CuppaPredictions cuppaPredictions)
    {
        Categories.ClfName mainCombinedClfName = cuppaPredictions.getMainCombinedClfName();

        CuppaPredictionEntry predictionEntry = cuppaPredictions
                .subsetByDataType(Categories.DataType.PROB)
                .subsetByClfName(mainCombinedClfName)
                .getTopPredictions(1)
                .get(0);

        return ImmutableProbabilityEntry.builder()
                .clfName(predictionEntry.ClfName.toString())
                .cancerType(predictionEntry.CancerType)
                .dataValue(predictionEntry.DataValue)
                .rank(predictionEntry.Rank)
                .rankGroup(predictionEntry.RankGroup)
                .build();
    }

    @NotNull
    public static List<ImmutableProbabilityEntry> getProbabilities(CuppaPredictions cuppaPredictions){

        CuppaPredictions probs = cuppaPredictions.subsetByDataType(Categories.DataType.PROB);

        List<ImmutableProbabilityEntry> probsNew = new ArrayList<>();
        for(CuppaPredictionEntry predictionEntry : probs.PredictionEntries)
        {
            ImmutableProbabilityEntry probabilityEntry = ImmutableProbabilityEntry.builder()
                    .clfName(predictionEntry.ClfName.toString())
                    .cancerType(predictionEntry.CancerType)
                    .dataValue(predictionEntry.DataValue)
                    .rank(predictionEntry.Rank)
                    .rankGroup(predictionEntry.RankGroup)
                    .build();
            probsNew.add(probabilityEntry);
        }
        return probsNew;
    }

    @NotNull
    public static List<ImmutableFeatureContributionEntry> getFeatureContributions(CuppaPredictions cuppaPredictions){

        CuppaPredictions contribs = cuppaPredictions.subsetByDataType(Categories.DataType.FEAT_CONTRIB);

        List<ImmutableFeatureContributionEntry> contribsNew = new ArrayList<>();
        for(CuppaPredictionEntry predictionEntry : contribs.PredictionEntries)
        {
            ImmutableFeatureContributionEntry contribsEntry = ImmutableFeatureContributionEntry.builder()
                    .clfName(predictionEntry.ClfName.toString())
                    .featName(predictionEntry.FeatName)
                    .featValue(predictionEntry.FeatValue)
                    .cancerType(predictionEntry.CancerType)
                    .dataValue(predictionEntry.DataValue)
                    .rank(predictionEntry.Rank)
                    .rankGroup(predictionEntry.RankGroup)
                    .build();

            contribsNew.add(contribsEntry);
        }
        return contribsNew;
    }

    @NotNull
    public static List<ImmutableSignatureQuantileEntry> getSignatureQuantiles(CuppaPredictions cuppaPredictions){

        CuppaPredictions sigQuantiles = cuppaPredictions.subsetByDataType(Categories.DataType.SIG_QUANTILE);

        List<ImmutableSignatureQuantileEntry> sigQuantilesNew = new ArrayList<>();
        for(CuppaPredictionEntry predictionEntry : sigQuantiles.PredictionEntries)
        {
            ImmutableSignatureQuantileEntry sigQuantilesEntry = ImmutableSignatureQuantileEntry.builder()
                    .featName(predictionEntry.FeatName)
                    .featValue(predictionEntry.FeatValue)
                    .cancerType(predictionEntry.CancerType)
                    .dataValue(predictionEntry.DataValue)
                    .rank(predictionEntry.Rank)
                    .rankGroup(predictionEntry.RankGroup)
                    .build();

            sigQuantilesNew.add(sigQuantilesEntry);
        }
        return sigQuantilesNew;
    }
}
