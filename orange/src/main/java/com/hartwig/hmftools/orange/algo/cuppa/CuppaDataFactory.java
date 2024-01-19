package com.hartwig.hmftools.orange.algo.cuppa;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

import com.hartwig.hmftools.common.cuppa2.Categories;
import com.hartwig.hmftools.common.cuppa2.CuppaPredictionEntry;
import com.hartwig.hmftools.common.cuppa2.CuppaPredictions;
import com.hartwig.hmftools.datamodel.cuppa.CuppaData;
import com.hartwig.hmftools.datamodel.cuppa.CuppaPrediction;
import com.hartwig.hmftools.datamodel.cuppa.ImmutableCuppaData;
import com.hartwig.hmftools.datamodel.cuppa.ImmutableCuppaPrediction;

import org.jetbrains.annotations.NotNull;

public final class CuppaDataFactory
{
    @NotNull
    public static CuppaData create(@NotNull CuppaPredictions cuppaPredictions) throws Exception
    {

        List<CuppaPrediction> predictions = extractSortedProbabilities(cuppaPredictions);

        CuppaPrediction bestPrediction = null;

        int NULL_INT_VALUE = -1;
        int lineCount = NULL_INT_VALUE;
        int maxComplexSize = NULL_INT_VALUE;
        int simpleDups32To200B = NULL_INT_VALUE;
        int telomericSGLs = NULL_INT_VALUE;

        if(predictions.size() != 0)
        {
            bestPrediction = predictions.get(0);
            lineCount = getSvFeatureValue(cuppaPredictions, "sv.LINE");
            maxComplexSize = getSvFeatureValue(cuppaPredictions, "sv.MAX_COMPLEX_SIZE");
            simpleDups32To200B = getSvFeatureValue(cuppaPredictions, "sv.SIMPLE_DEL_20KB_1MB");
            telomericSGLs = getSvFeatureValue(cuppaPredictions, "sv.TELOMERIC_SGL");
        }

        return ImmutableCuppaData.builder()
                .predictions(predictions)
                .bestPrediction(bestPrediction)
                .lineCount(lineCount)
                .maxComplexSize(maxComplexSize)
                .simpleDups32To200B(simpleDups32To200B)
                .telomericSGLs(telomericSGLs)
                .build();
    }

    @NotNull
    public static List<CuppaPrediction> extractSortedProbabilities(@NotNull CuppaPredictions cuppaPredictions)
    {
        CuppaPredictions probabilitiesAllClassifiers = cuppaPredictions
                .subsetByDataType(Categories.DataType.PROB)
                .sortByRank();

        List<CuppaPrediction> cuppaPredictionsOrangeFormat = new ArrayList<>();
        for(String cancerType : probabilitiesAllClassifiers.CancerTypes)
        {
            CuppaPredictions probabilitiesOneCancerType = probabilitiesAllClassifiers.subsetByCancerType(cancerType);

            Map<Categories.ClfName, Double> probabilitiesByClassifier = new HashMap<>();
            for(int i = 0; i < probabilitiesOneCancerType.size(); i++)
            {
                probabilitiesByClassifier.put(
                        probabilitiesOneCancerType.get(i).ClfName,
                        probabilitiesOneCancerType.get(i).DataValue
                );
            }

            CuppaPrediction prediction = ImmutableCuppaPrediction.builder()
                    .cancerType(cancerType)
                    .likelihood(probabilitiesByClassifier.get(probabilitiesAllClassifiers.MainCombinedClfName))
                    .genomicPositionClassifier(probabilitiesByClassifier.get(Categories.ClfName.GEN_POS))
                    .snvPairwiseClassifier(probabilitiesByClassifier.get(Categories.ClfName.SNV96))
                    .featureClassifier(probabilitiesByClassifier.get(Categories.ClfName.EVENT))
                    .expressionPairwiseClassifier(probabilitiesByClassifier.get(Categories.ClfName.GENE_EXP))
                    .altSjCohortClassifier(probabilitiesByClassifier.get(Categories.ClfName.ALT_SJ))
                    .build();

            cuppaPredictionsOrangeFormat.add(prediction);
        }
        return cuppaPredictionsOrangeFormat;
    }

    public static int getSvFeatureValue(CuppaPredictions cuppaPredictions, String featureName) throws Exception
    {
        // Feature values are replicated for each cancer type because the `cuppaPredictions` table is in long form.
        // Use `.findFirst()` to get a single value
        CuppaPredictionEntry predictionEntry = cuppaPredictions.PredictionEntries.stream()
                .filter(o -> o.DataType == Categories.DataType.FEAT_CONTRIB & o.FeatName.equals(featureName))
                .findFirst()
                .orElseThrow(() -> new Exception("Input CuppaPredictions is empty"));

        return (int) Math.round(predictionEntry.FeatValue);
    }
}
