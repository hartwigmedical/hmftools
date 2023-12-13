package com.hartwig.hmftools.orange.algo.cuppa2;

import java.util.ArrayList;
import java.util.List;

import com.hartwig.hmftools.common.cuppa2.Categories;
import com.hartwig.hmftools.common.cuppa2.CuppaVisDataEntry;
import com.hartwig.hmftools.common.cuppa2.CuppaVisData;
import com.hartwig.hmftools.datamodel.cuppa2.ImmutableCuppa2Data;
import com.hartwig.hmftools.datamodel.cuppa2.ImmutableProbabilityEntry;
import com.hartwig.hmftools.datamodel.cuppa2.ImmutableFeatureContributionEntry;
import com.hartwig.hmftools.datamodel.cuppa2.ImmutableSignatureQuantileEntry;

import org.jetbrains.annotations.NotNull;

public class Cuppa2DataFactory
{
    @NotNull
    public static ImmutableCuppa2Data create(@NotNull CuppaVisData visData)
    {
        return ImmutableCuppa2Data.builder()
                .topPrediction(getTopPrediction(visData))
                .probs(getProbabilities(visData))
                .featContribs(getFeatureContributions(visData))
                .sigQuantiles(getSignatureQuantiles(visData))
                .build();
    }

    @NotNull
    public static ImmutableProbabilityEntry getTopPrediction(CuppaVisData visData)
    {
        Categories.ClfName mainCombinedClfName = visData.getMainCombinedClfName();

        CuppaVisDataEntry visDataEntry = visData
                .subsetByDataType(Categories.DataType.PROB)
                .subsetByClfName(mainCombinedClfName)
                .getTopPredictions(1)
                .get(0);

        return ImmutableProbabilityEntry.builder()
                .clfName(visDataEntry.ClfName.toString())
                .cancerType(visDataEntry.CancerType)
                .dataValue(visDataEntry.DataValue)
                .rank(visDataEntry.Rank)
                .rankGroup(visDataEntry.RankGroup)
                .build();
    }

    @NotNull
    public static List<ImmutableProbabilityEntry> getProbabilities(CuppaVisData visData){

        CuppaVisData probs = visData.subsetByDataType(Categories.DataType.PROB);

        List<ImmutableProbabilityEntry> probsNew = new ArrayList<>();
        for(CuppaVisDataEntry visDataEntry : probs.VisDataEntries)
        {
            ImmutableProbabilityEntry probabilityEntry = ImmutableProbabilityEntry.builder()
                    .clfName(visDataEntry.ClfName.toString())
                    .cancerType(visDataEntry.CancerType)
                    .dataValue(visDataEntry.DataValue)
                    .rank(visDataEntry.Rank)
                    .rankGroup(visDataEntry.RankGroup)
                    .build();
            probsNew.add(probabilityEntry);
        }
        return probsNew;
    }

    @NotNull
    public static List<ImmutableFeatureContributionEntry> getFeatureContributions(CuppaVisData visData){

        CuppaVisData contribs = visData.subsetByDataType(Categories.DataType.FEAT_CONTRIB);

        List<ImmutableFeatureContributionEntry> contribsNew = new ArrayList<>();
        for(CuppaVisDataEntry visDataEntry : contribs.VisDataEntries)
        {
            ImmutableFeatureContributionEntry contribsEntry = ImmutableFeatureContributionEntry.builder()
                    .clfName(visDataEntry.ClfName.toString())
                    .featName(visDataEntry.FeatName)
                    .featValue(visDataEntry.FeatValue)
                    .cancerType(visDataEntry.CancerType)
                    .dataValue(visDataEntry.DataValue)
                    .rank(visDataEntry.Rank)
                    .rankGroup(visDataEntry.RankGroup)
                    .build();

            contribsNew.add(contribsEntry);
        }
        return contribsNew;
    }

    @NotNull
    public static List<ImmutableSignatureQuantileEntry> getSignatureQuantiles(CuppaVisData visData){

        CuppaVisData sigQuantiles = visData.subsetByDataType(Categories.DataType.SIG_QUANTILE);

        List<ImmutableSignatureQuantileEntry> sigQuantilesNew = new ArrayList<>();
        for(CuppaVisDataEntry visDataEntry : sigQuantiles.VisDataEntries)
        {
            ImmutableSignatureQuantileEntry sigQuantilesEntry = ImmutableSignatureQuantileEntry.builder()
                    .featName(visDataEntry.FeatName)
                    .featValue(visDataEntry.FeatValue)
                    .cancerType(visDataEntry.CancerType)
                    .dataValue(visDataEntry.DataValue)
                    .rank(visDataEntry.Rank)
                    .rankGroup(visDataEntry.RankGroup)
                    .build();

            sigQuantilesNew.add(sigQuantilesEntry);
        }
        return sigQuantilesNew;
    }
}
