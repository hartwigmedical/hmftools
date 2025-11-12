package com.hartwig.hmftools.qsee.cohort;

import static com.hartwig.hmftools.qsee.common.QseeConstants.QC_LOGGER;

import static org.junit.Assert.assertArrayEquals;
import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertTrue;

import java.util.Arrays;

import org.apache.logging.log4j.Level;
import org.apache.logging.log4j.core.config.Configurator;
import org.junit.Before;
import org.junit.Test;

public class PercentileTransformerTest
{
    private final int DEFAULT_PERCENTILE_INCREMENT = 10;
    private final double[] DEFAULT_COHORT_VALUES = { 0, 0, 5, 10, 10, 10, 10, 10, 10, 10, Double.NaN, Double.POSITIVE_INFINITY };

    @Before
    public void setUp()
    {
        Configurator.setLevel(QC_LOGGER.getName(), Level.ERROR);
    }

    private PercentileTransformer createDefaultTransformer()
    {
        PercentileTransformer transformer = PercentileTransformer.withIncrement(DEFAULT_PERCENTILE_INCREMENT);
        transformer.fit(DEFAULT_COHORT_VALUES);
        return transformer;
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
    public void nanAndInfIgnoredDuringFit()
    {
        PercentileTransformer transformerWithNan = createDefaultTransformer();

        PercentileTransformer transformerWithoutNan = PercentileTransformer.withIncrement(DEFAULT_PERCENTILE_INCREMENT);
        double[] cohortValuesWithoutNan = Arrays.stream(DEFAULT_COHORT_VALUES).filter(value -> !Double.isNaN(value)).toArray();
        transformerWithoutNan.fit(cohortValuesWithoutNan);

        assertArrayEquals(transformerWithNan.getPercentiles(), transformerWithoutNan.getPercentiles(), 0.001);
        assertArrayEquals(transformerWithNan.getRefValues(), transformerWithoutNan.getRefValues(), 0.001);
    }

    @Test
    public void cohortValuesDedupedDuringFit()
    {
        PercentileTransformer transformer = createDefaultTransformer();

        double[] expectedDedupPercentiles = { 5, 20, 30, 70 };
        double[] expectedDedupRefValues = { 0, 4, 8.5, 10 };

        double[] actualDedupPercentiles = transformer.getPercentilesDeduped();
        double[] actualDedupRefValues = transformer.getRefValuesDeduped();

        assertArrayEquals(expectedDedupPercentiles, actualDedupPercentiles, 0.001);
        assertArrayEquals(expectedDedupRefValues, actualDedupRefValues, 0.001);
    }

    @Test
    public void canInterpolatePercentileToFeatureValue()
    {
        PercentileTransformer transformer = createDefaultTransformer();
        assertEquals(0, transformer.percentileToFeatureValue(0), 0.001);
        assertEquals((4 + 8.5) / 2, transformer.percentileToFeatureValue(25), 0.001);
    }

    @Test
    public void canInterpolateFeatureValuesToPercentiles()
    {
        PercentileTransformer transformer = createDefaultTransformer();

        // Matches ref data exactly
        assertEquals(5, transformer.featureValueToPercentile(0), 0.01);
        assertEquals(20, transformer.featureValueToPercentile(4), 0.01);
        assertEquals(30, transformer.featureValueToPercentile(8.5), 0.01);
        assertEquals(70, transformer.featureValueToPercentile(10), 0.01);

        // Requires interpolation
        assertEquals(12.50, transformer.featureValueToPercentile(2), 0.01);
        assertEquals(26.67, transformer.featureValueToPercentile(7), 0.01);

        // Out of bounds
        assertEquals(Double.POSITIVE_INFINITY, transformer.featureValueToPercentile(1000), 0.01);
        assertEquals(Double.NEGATIVE_INFINITY, transformer.featureValueToPercentile(-1000), 0.01);
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
        double actualPercentile = transformer.featureValueToPercentile(0);

        assertEquals(expectedPercentile, actualPercentile, 0.1);
    }

    @Test
    public void canCreateFromPrefitData()
    {
        double[] prefitPercentiles = { 0, 25, 50, 75, 100 };
        double[] prefitRefValues = { 0, 0, 10, 10, 10 };

        PercentileTransformer newTransformer = PercentileTransformer.fromPrefitData(prefitPercentiles, prefitRefValues);

        double[] expectedDedupPercentiles = { 12.5, 75 };
        double[] expectedDedupRefValues = { 0, 10 };

        double[] actualDedupPercentiles = newTransformer.getPercentilesDeduped();
        double[] actualDedupRefValues = newTransformer.getRefValuesDeduped();

        assertArrayEquals(expectedDedupPercentiles, actualDedupPercentiles, 0.01);
        assertArrayEquals(expectedDedupRefValues, actualDedupRefValues, 0.01);
    }
}
