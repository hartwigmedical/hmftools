package com.hartwig.hmftools.finding.util;

import java.util.List;
import java.util.Map;

import com.hartwig.hmftools.finding.FindingUtil;
import com.hartwig.hmftools.finding.datamodel.FindingRecord;
import com.hartwig.hmftools.finding.datamodel.FindingRecordBuilder;
import com.hartwig.hmftools.finding.datamodel.PredictedTumorOrigin;
import com.hartwig.hmftools.finding.datamodel.PredictedTumorOriginBuilder;
import com.hartwig.hmftools.finding.datamodel.PredictedTumorOriginPredictionBuilder;
import com.hartwig.hmftools.finding.datamodel.finding.FindingItem;
import com.hartwig.hmftools.finding.datamodel.finding.FindingItemBuilder;
import com.hartwig.hmftools.finding.datamodel.finding.FindingStatus;

import jakarta.validation.constraints.NotNull;

public class PTOConverter
{
    // TODO: Is this constant defined elsewhere?
    private static final double CUPPA_INCONCLUSIVE_CUT_OFF = 0.8;
    private static final double BEST_LIKELIHOOD_CUT_OFF = 0.5;

    private static final String LOWER_GI_TRACT = "Lower GI tract";
    // TODO: Check if these mappings are still necessary
    private static final Map<String, String> CURATED_CANCER_TYPES =
            Map.of("Uterus: Endometrium", "Endometrium", "Colorectum/Appendix/SmallIntestine", LOWER_GI_TRACT, "Colorectum/Small intestine/Appendix", LOWER_GI_TRACT);

    @NotNull
    public static FindingRecord convert(@NotNull FindingRecord record)
    {
        return FindingRecordBuilder.builder(record)
                .predictedTumorOrigin(convert(record.predictedTumorOrigin()))
                .build();
    }

    @NotNull
    private static FindingItem<PredictedTumorOrigin> convert(@NotNull FindingItem<PredictedTumorOrigin> findingItem)
    {
        PredictedTumorOrigin predictedTumorOrigin = findingItem.finding();
        if(predictedTumorOrigin != null)
        {
            List<PredictedTumorOrigin.Prediction> predictions = predictedTumorOrigin.predictions();
            if(!predictions.isEmpty())
            {
                PredictedTumorOrigin.Prediction best = predictions.get(0);
                predictions = predictions.stream()
                        .filter(prediction ->
                                prediction.likelihood() >= CUPPA_INCONCLUSIVE_CUT_OFF)
                        .map(p -> PredictedTumorOriginPredictionBuilder.builder(p).cancerType(curateCancerType(p.cancerType()))
                                .build())
                        .toList();
                FindingStatus findingStatus = findingItem.status();
                if(predictions.isEmpty())
                {
                    // Changing status code because this is different from there being no results.
                    // The issue is that no results meet the required criteria.

                    if (best.likelihood() >= BEST_LIKELIHOOD_CUT_OFF)
                    {
                        findingStatus = FindingUtil.noReportableValueStatus(findingStatus, FindingStatus.Status.NOT_RELIABLE);
                        predictions = List.of(best);
                    } else {
                        findingStatus = FindingUtil.noReportableValueStatus(findingStatus, FindingStatus.Status.NOT_AVAILABLE);
                    }
                }
                return FindingItemBuilder.builder(findingItem)
                        .status(findingStatus)
                        .finding(PredictedTumorOriginBuilder.builder(predictedTumorOrigin)
                                .predictions(predictions)
                                .build())
                        .build();
            }
        }
        return findingItem;
    }

    @NotNull
    private static String curateCancerType(@NotNull String cancerType)
    {
        return CURATED_CANCER_TYPES.getOrDefault(cancerType, cancerType);
    }
}
