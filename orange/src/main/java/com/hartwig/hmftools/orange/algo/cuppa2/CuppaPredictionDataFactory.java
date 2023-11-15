package com.hartwig.hmftools.orange.algo.cuppa2;

import java.util.ArrayList;
import java.util.List;

import com.hartwig.hmftools.common.cuppa.CuppaDataFile;
import com.hartwig.hmftools.common.cuppa.SvDataType;
import com.hartwig.hmftools.common.cuppa2.Categories;
import com.hartwig.hmftools.common.cuppa2.CuppaPredictionEntry;
import com.hartwig.hmftools.common.cuppa2.CuppaPredictions;
import com.hartwig.hmftools.datamodel.cuppa2.ImmutableCuppaPredictionData;
import com.hartwig.hmftools.datamodel.cuppa2.CuppaPredictionData;

import org.jetbrains.annotations.NotNull;

public class CuppaPredictionDataFactory
{
    public static CuppaPredictionData create(@NotNull List<CuppaDataFile> entries)
    {
        return ImmutableCuppaPredictionData.builder()
                .predictions(extractPredictions(entries))
                .simpleDups32To200B(safeInt(entries, SvDataType.SIMPLE_DUP_32B_200B))
                .maxComplexSize(safeInt(entries, SvDataType.MAX_COMPLEX_SIZE))
                .lineCount(safeInt(entries, SvDataType.LINE))
                .telomericSGLs(safeInt(entries, SvDataType.TELOMERIC_SGL))
                .build();
    }

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

    public static List<ProbabilityEntry> getProbabilities(CuppaPredictions cuppaPredictions){

        CuppaPredictions probs = cuppaPredictions.subsetByDataType(Categories.DataType.PROB);

        List<ProbabilityEntry> probsNew = new ArrayList<>();
        for(CuppaPredictionEntry predictionEntry : probs.PredictionEntries)
        {
            probsNew.add(ProbabilityEntry.fromCuppaPredictionEntry(predictionEntry));
        }
        return probsNew;
    }

    public static List<FeatureContributionEntry> getFeatureContributions(CuppaPredictions cuppaPredictions){

        CuppaPredictions contribs = cuppaPredictions.subsetByDataType(Categories.DataType.FEAT_CONTRIB);

        List<FeatureContributionEntry> contribsNew = new ArrayList<>();
        for(CuppaPredictionEntry predictionEntry : contribs.PredictionEntries)
        {
            contribsNew.add(FeatureContributionEntry.fromCuppaPredictionEntry(predictionEntry));
        }
        return contribsNew;
    }

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
