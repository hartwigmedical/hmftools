package com.hartwig.hmftools.finding.util;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertFalse;
import static org.junit.Assert.assertNotNull;
import static org.junit.Assert.assertNull;
import static org.junit.Assert.assertTrue;

import java.util.List;
import java.util.Set;

import com.hartwig.hmftools.finding.datamodel.FindingRecord;
import com.hartwig.hmftools.finding.datamodel.PredictedTumorOrigin;
import com.hartwig.hmftools.finding.datamodel.TestFindingFactory;
import com.hartwig.hmftools.finding.datamodel.TestFindingRecordFactory;
import com.hartwig.hmftools.finding.datamodel.finding.FindingItem;
import com.hartwig.hmftools.finding.datamodel.finding.FindingStatus;

import org.junit.Test;

public class PTOConverterTest
{
    @Test
    public void convertWithResults()
    {
        FindingRecord original = createPTOFindingItem(0.8);
        FindingRecord converted = PTOConverter.convert(original);

        FindingItem<PredictedTumorOrigin> findingItem = converted.predictedTumorOrigin();
        assertEquals(FindingStatus.Status.OK, findingItem.status().status());
        assertEquals(Set.of(), findingItem.status().errors());

        PredictedTumorOrigin predictedTumorOrigin = findingItem.finding();
        assertNotNull(predictedTumorOrigin);
        assertNotNull(predictedTumorOrigin.bestPredictionLikelihood());
        assertEquals(0.8, predictedTumorOrigin.bestPredictionLikelihood(), 0.00001);
        assertFalse(predictedTumorOrigin.predictions().isEmpty());
    }

    @Test
    public void convertNoResultsWithBestLikelihood()
    {
        FindingRecord original = createPTOFindingItem(0.7);
        FindingRecord converted = PTOConverter.convert(original);

        FindingItem<PredictedTumorOrigin> findingItem = converted.predictedTumorOrigin();
        assertEquals(FindingStatus.Status.NOT_AVAILABLE, findingItem.status().status());
        assertEquals(Set.of(FindingStatus.Issue.NO_REPORTABLE_VALUE), findingItem.status().errors());

        PredictedTumorOrigin predictedTumorOrigin = findingItem.finding();
        assertNotNull(predictedTumorOrigin);
        assertNotNull(predictedTumorOrigin.bestPredictionLikelihood());
        assertEquals(0.7, predictedTumorOrigin.bestPredictionLikelihood(), 0.00001);
        assertTrue(predictedTumorOrigin.predictions().isEmpty());
    }

    @Test
    public void convertNoResultsWithoutBestLikelihood()
    {
        FindingRecord original = createPTOFindingItem(0.4);
        FindingRecord converted = PTOConverter.convert(original);

        FindingItem<PredictedTumorOrigin> findingItem = converted.predictedTumorOrigin();
        assertEquals(FindingStatus.Status.NOT_AVAILABLE, findingItem.status().status());
        assertEquals(Set.of(FindingStatus.Issue.NO_REPORTABLE_VALUE), findingItem.status().errors());

        PredictedTumorOrigin predictedTumorOrigin = findingItem.finding();
        assertNotNull(predictedTumorOrigin);
        assertNull(predictedTumorOrigin.bestPredictionLikelihood());
        assertTrue(predictedTumorOrigin.predictions().isEmpty());
    }

    private static FindingRecord createPTOFindingItem(double likelihood)
    {
        return TestFindingRecordFactory.createMinimalTestFindingRecordBuilder()
                .predictedTumorOrigin(TestFindingFactory.buildFindingItem(FindingStatus.Status.OK, TestFindingFactory.predictedTumorOriginBuilder()
                        .predictions(List.of(TestFindingFactory.predictedTumorOriginPredictionBuilder()
                                .likelihood(likelihood)
                                .build()))
                        .bestPredictionLikelihood(likelihood)
                        .build()))
                .build();
    }
}
