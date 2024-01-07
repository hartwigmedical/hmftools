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
    public static CuppaData create(@NotNull CuppaPredictions cuppaPredictions){

        List<CuppaPrediction> predictions = extractSortedProbabilities(cuppaPredictions);

        CuppaPrediction bestPrediction = null;
        Integer lineCount = null;
        Integer maxComplexSize = null;
        Integer simpleDups32To200B = null;
        Integer telomericSGLs = null;

        if(predictions.size() != 0)
        {
            bestPrediction = predictions.get(0);
            lineCount = getFeatureValue(cuppaPredictions, "sv.LINE").intValue();
            maxComplexSize = getFeatureValue(cuppaPredictions, "sv.MAX_COMPLEX_SIZE").intValue();
            simpleDups32To200B = getFeatureValue(cuppaPredictions, "sv.SIMPLE_DEL_20KB_1MB").intValue();
            telomericSGLs = getFeatureValue(cuppaPredictions, "sv.TELOMERIC_SGL").intValue();
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
                    .cuppaMajorVersion("v2")
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

    @NotNull
    public static Double getFeatureValue(CuppaPredictions cuppaPredictions, String featureName)
    {
        // Feature values are replicated for each cancer type because the `cuppaPredictions` table is in long form.
        // Use `.findFirst()` to get a single value
        CuppaPredictionEntry predictionEntry = cuppaPredictions.PredictionEntries.stream()
                .filter(o -> o.DataType == Categories.DataType.FEAT_CONTRIB & o.FeatName.equals(featureName))
                .findFirst()
                .get();

        return predictionEntry.FeatValue;
    }
}
