package com.hartwig.hmftools.qsee.lod;

import static org.junit.Assert.assertEquals;

import java.util.Map;

import com.hartwig.hmftools.common.metrics.DepthRange;

import org.junit.Test;

public class LimitOfDetectionCalcsTest
{
    @Test
    public void canCalcTumorVafLod()
    {
        double actualVaf;
        double expectedVaf;

        Map<DepthRange,Integer> depthCounts = Map.of(
                new DepthRange(20, 20), 10,
                new DepthRange(100, 100), 90
        );

        // Without purity adjustment
        actualVaf = LimitOfDetectionCalcs.calcMinTumorVafToDetectVariant(depthCounts, 1.0);
        expectedVaf = 0.245;
        assertEquals(expectedVaf, actualVaf, 1e-3);

        // With purity adjustment
        actualVaf = LimitOfDetectionCalcs.calcMinTumorVafToDetectVariant(depthCounts, 0.5);
        expectedVaf = 0.490;
        assertEquals(expectedVaf, actualVaf, 1e-3);

        // Purity adjustment leads to VAF > 1 (not achievable)
        actualVaf = LimitOfDetectionCalcs.calcMinTumorVafToDetectVariant(depthCounts, 0.2);
        expectedVaf = 1.224;
        assertEquals(expectedVaf, actualVaf, 1e-3);

        // Depth too low, no VAF can achieve the missed variant likelihood target
        Map<DepthRange,Integer> lowDepthCounts = Map.of(
                new DepthRange(2, 2), 10,
                new DepthRange(100, 100), 90
        );

        actualVaf = LimitOfDetectionCalcs.calcMinTumorVafToDetectVariant(lowDepthCounts, 1.0);
        expectedVaf = 1.0;
        assertEquals(expectedVaf, actualVaf, 0);
    }
}
