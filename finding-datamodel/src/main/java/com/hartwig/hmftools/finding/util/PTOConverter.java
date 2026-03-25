package com.hartwig.hmftools.finding.util;

import static com.hartwig.hmftools.finding.datamodel.finding.FindingStatus.Issue.NO_REPORTABLE_VALUE;

import java.util.List;
import java.util.Set;
import java.util.TreeSet;

import com.hartwig.hmftools.finding.datamodel.FindingRecord;
import com.hartwig.hmftools.finding.datamodel.FindingRecordBuilder;
import com.hartwig.hmftools.finding.datamodel.PredictedTumorOrigin;
import com.hartwig.hmftools.finding.datamodel.PredictedTumorOriginBuilder;
import com.hartwig.hmftools.finding.datamodel.finding.FindingItem;
import com.hartwig.hmftools.finding.datamodel.finding.FindingItemBuilder;
import com.hartwig.hmftools.finding.datamodel.finding.FindingStatus;
import com.hartwig.hmftools.finding.datamodel.finding.FindingStatusBuilder;

import jakarta.validation.constraints.NotNull;

public class PTOConverter
{
    // TODO: Is this constant defined elsewhere?
    private static final double CUPPA_INCONCLUSIVE_CUT_OFF = 0.8;

    public static FindingRecord convert(FindingRecord record)
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
                predictions = predictions.stream().filter(prediction ->
                        prediction.likelihood() >= CUPPA_INCONCLUSIVE_CUT_OFF).toList();
                if(predictions.isEmpty())
                {
                    // Changing status code because this is different from there being no results.
                    // The issue is that no results meet the required criteria.
                    FindingStatusBuilder.builder(findingItem.status())
                            .status(FindingStatus.Status.NOT_AVAILABLE)
                            .errors(new TreeSet<>(Set.of(NO_REPORTABLE_VALUE)))
                            .build();
                }
                return FindingItemBuilder.builder(findingItem)
                        .finding(PredictedTumorOriginBuilder.builder(predictedTumorOrigin)
                                .predictions(predictions)
                                .build())
                        .build();
            }
        }
        return findingItem;
    }
}
