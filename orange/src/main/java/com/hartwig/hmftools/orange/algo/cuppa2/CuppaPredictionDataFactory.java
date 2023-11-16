package com.hartwig.hmftools.orange.algo.cuppa2;

import java.util.ArrayList;
import java.util.List;

import com.hartwig.hmftools.common.cuppa2.Categories;
import com.hartwig.hmftools.common.cuppa2.CuppaPredictionEntry;
import com.hartwig.hmftools.common.cuppa2.CuppaPredictions;
import com.hartwig.hmftools.common.cuppa2.FeatureContributionEntry;
import com.hartwig.hmftools.common.cuppa2.ProbabilityEntry;
import com.hartwig.hmftools.common.cuppa2.SignatureQuantileEntry;
import com.hartwig.hmftools.datamodel.cuppa2.ImmutableCuppaPredictionData;

import org.jetbrains.annotations.NotNull;

public class CuppaPredictionDataFactory
{
    @NotNull
    public static ImmutableCuppaPredictionData create(@NotNull CuppaPredictions cuppaPredictions)
    {
        return ImmutableCuppaPredictionData.builder()
                .topPrediction(getTopPrediction(cuppaPredictions))
                .probs(getProbabilities(cuppaPredictions))
                .featContribs(getFeatureContributions(cuppaPredictions))
                .sigQuantiles(getSignatureQuantiles(cuppaPredictions))
                .build();
    }

    @NotNull
    public static ProbabilityEntry getTopPrediction(CuppaPredictions cuppaPredictions)
    {
        Categories.ClfName mainCombinedClfName = cuppaPredictions.getMainCombinedClfName();

        CuppaPredictionEntry predictionEntry = cuppaPredictions
                .subsetByDataType(Categories.DataType.PROB)
                .subsetByClfName(mainCombinedClfName)
                .getTopPredictions(1)
                .get(0);

        return ProbabilityEntry.fromCuppaPredictionEntry(predictionEntry);
    }

    @NotNull
    public static List<ProbabilityEntry> getProbabilities(CuppaPredictions cuppaPredictions){

        CuppaPredictions probs = cuppaPredictions.subsetByDataType(Categories.DataType.PROB);

        List<ProbabilityEntry> probsNew = new ArrayList<>();
        for(CuppaPredictionEntry predictionEntry : probs.PredictionEntries)
        {
            probsNew.add(ProbabilityEntry.fromCuppaPredictionEntry(predictionEntry));
        }
        return probsNew;
    }

    @NotNull
    public static List<FeatureContributionEntry> getFeatureContributions(CuppaPredictions cuppaPredictions){

        CuppaPredictions contribs = cuppaPredictions.subsetByDataType(Categories.DataType.FEAT_CONTRIB);

        List<FeatureContributionEntry> contribsNew = new ArrayList<>();
        for(CuppaPredictionEntry predictionEntry : contribs.PredictionEntries)
        {
            contribsNew.add(FeatureContributionEntry.fromCuppaPredictionEntry(predictionEntry));
        }
        return contribsNew;
    }

    @NotNull
    public static List<SignatureQuantileEntry> getSignatureQuantiles(CuppaPredictions cuppaPredictions){

        CuppaPredictions sigQuantiles = cuppaPredictions.subsetByDataType(Categories.DataType.SIG_QUANTILE);

        List<SignatureQuantileEntry> sigQuantilesNew = new ArrayList<>();
        for(CuppaPredictionEntry predictionEntry : sigQuantiles.PredictionEntries)
        {
            sigQuantilesNew.add(SignatureQuantileEntry.fromCuppaPredictionEntry(predictionEntry));
        }
        return sigQuantilesNew;
    }
}
