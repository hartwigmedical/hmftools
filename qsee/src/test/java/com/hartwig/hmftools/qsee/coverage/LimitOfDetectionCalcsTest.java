package com.hartwig.hmftools.qsee.coverage;

import static org.junit.Assert.assertEquals;

import java.util.Map;

import com.hartwig.hmftools.common.metrics.DepthRange;

import org.junit.Test;

public class LimitOfDetectionCalcsTest
{
    @Test
    public void canCalcLodVaf()
    {
        double actualVaf;
        double expectedVaf;

        Map<DepthRange,Integer> depthCounts = Map.of(
                new DepthRange(20, 20), 10,
                new DepthRange(100, 100), 90
        );

        actualVaf = LimitOfDetectionCalcs.calcMinSampleVafToDetectVariant(depthCounts);
        expectedVaf = 0.132;
        assertEquals(expectedVaf, actualVaf, 1e-3);

        // Depth too low, no VAF can achieve the missed variant likelihood target
        Map<DepthRange,Integer> lowDepthCounts = Map.of(
                new DepthRange(2, 2), 10,
                new DepthRange(100, 100), 90
        );

        actualVaf = LimitOfDetectionCalcs.calcMinSampleVafToDetectVariant(lowDepthCounts);
        expectedVaf = 1.0;
        assertEquals(expectedVaf, actualVaf, 0);

        // Depth high enough to give LoD VAF below lower bound
        Map<DepthRange,Integer> highDepthCounts = Map.of(
                new DepthRange(1000, 1000), 100
        );

        actualVaf = LimitOfDetectionCalcs.calcMinSampleVafToDetectVariant(highDepthCounts);
        expectedVaf = LimitOfDetectionCalcs.LOD_VAF_LOWER_BOUND;
        assertEquals(expectedVaf, actualVaf, 0);
    }
}
