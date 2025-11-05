package com.hartwig.hmftools.qsee.cohort;

import static org.junit.Assert.assertArrayEquals;
import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertTrue;

import java.util.Arrays;

import org.junit.Test;

public class PercentileTransformerTest
{
    private final int DEFAULT_PERCENTILE_INCREMENT = 10;
    private final double[] DEFAULT_COHORT_VALUES = { 0, 0, 5, 10, 10, 10, 10, 10, 10, 10, Double.NaN, Double.POSITIVE_INFINITY };

    @Test
    public void nanAndInfIgnoredDuringFit()
    {
        PercentileTransformer transformerWithNan = PercentileTransformer.withIncrement(DEFAULT_PERCENTILE_INCREMENT);
        transformerWithNan.fit(DEFAULT_COHORT_VALUES);

        PercentileTransformer transformerWithoutNan = PercentileTransformer.withIncrement(DEFAULT_PERCENTILE_INCREMENT);
        double[] cohortValuesWithoutNan = Arrays.stream(DEFAULT_COHORT_VALUES).filter(value -> !Double.isNaN(value)).toArray();
        transformerWithoutNan.fit(cohortValuesWithoutNan);

        assertArrayEquals(transformerWithNan.getPercentiles(), transformerWithoutNan.getPercentiles(), 0.001);
        assertArrayEquals(transformerWithNan.getRefValues(), transformerWithoutNan.getRefValues(), 0.001);
    }

    @Test
    public void cohortValuesDedupedDuringFit()
    {
        PercentileTransformer transformer = PercentileTransformer.withIncrement(DEFAULT_PERCENTILE_INCREMENT);
        transformer.fit(DEFAULT_COHORT_VALUES);

        double[] expectedDedupPercentiles = { 5, 20, 30, 70 };
        double[] expectedDedupRefValues = { 0, 4, 8.5, 10 };

        double[] actualDedupPercentiles = transformer.getPercentilesDeduped();
        double[] actualDedupRefValues = transformer.getRefValuesDeduped();

        assertArrayEquals(expectedDedupPercentiles, actualDedupPercentiles, 0.001);
        assertArrayEquals(expectedDedupRefValues, actualDedupRefValues, 0.001);
    }

    @Test
    public void canTransformCohortValuesToPercentiles()
    {
        PercentileTransformer transformer = PercentileTransformer.withIncrement(DEFAULT_PERCENTILE_INCREMENT);
        transformer.fit(DEFAULT_COHORT_VALUES);

        // Matches ref data exactly
        assertEquals(5,  transformer.transform(0),   0.01);
        assertEquals(20, transformer.transform(4),   0.01);
        assertEquals(30, transformer.transform(8.5), 0.01);
        assertEquals(70, transformer.transform(10),  0.01);

        // Requires interpolation
        assertEquals(12.50, transformer.transform(2),  0.01);
        assertEquals(26.67, transformer.transform(7),  0.01);

        // Out of bounds
        assertEquals(Double.POSITIVE_INFINITY, transformer.transform( 1000),  0.01);
        assertEquals(Double.NEGATIVE_INFINITY, transformer.transform(-1000),  0.01);
    }

    @Test
    public void canFitOneCohortValue()
    {
        PercentileTransformer transformer = PercentileTransformer.withNumPercentiles(3);
        double[] cohortValues = { 5 };
        transformer.fit(cohortValues);

        double[] expectedPercentiles = { 50 };
        double[] expectedRefValues = { 5 };

        double[] actualPercentiles = transformer.getPercentilesDeduped();
        double[] actualRefValues = transformer.getRefValuesDeduped();

        assertArrayEquals(expectedRefValues, actualRefValues, 0.1);
        assertArrayEquals(expectedPercentiles, actualPercentiles, 0.1);
    }

    @Test
    public void canFitAndTransformFromAllNanCohortValues()
    {
        PercentileTransformer transformer = PercentileTransformer.withNumPercentiles(3);
        double[] cohortValues = { Double.NaN, Double.NaN };

        transformer.fit(cohortValues);

        double expectedPercentile = Double.NaN;
        double actualPercentile = transformer.transform(0);

        assertEquals(expectedPercentile, actualPercentile, 0.1);
    }

    @Test
    public void canCreateFromNumPercentilesOrIncrement()
    {
        double[] expectedPercentiles = { 0, 25, 50, 75, 100 };

        PercentileTransformer transformerNumPct = PercentileTransformer.withNumPercentiles(5);
        PercentileTransformer transformerIncrement = PercentileTransformer.withIncrement(25);

        assertArrayEquals(expectedPercentiles, transformerNumPct.getPercentiles(), 0.1);
        assertArrayEquals(expectedPercentiles, transformerIncrement.getPercentiles(), 0.1);
    }

    @Test
    public void canCreateFromPrefitData()
    {
        double[] prefitPercentiles = { 0, 25, 50, 75, 100 };
        double[] prefitRefValues = { 0, 0, 10, 10, 10 };

        PercentileTransformer newTransformer = PercentileTransformer.fromPrefitData(prefitPercentiles, prefitRefValues);

        double[] expectedDedupPercentiles = { 12.5, 75 };
        double[] expectedDedupRefValues   = {    0, 10 };

        double[] actualDedupPercentiles = newTransformer.getPercentilesDeduped();
        double[] actualDedupRefValues = newTransformer.getRefValuesDeduped();

        assertArrayEquals(expectedDedupPercentiles, actualDedupPercentiles, 0.01);
        assertArrayEquals(expectedDedupRefValues, actualDedupRefValues, 0.01);
    }
}
