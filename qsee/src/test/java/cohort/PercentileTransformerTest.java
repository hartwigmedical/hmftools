package cohort;

import static org.junit.Assert.assertArrayEquals;
import static org.junit.Assert.assertEquals;

import org.junit.Before;
import org.junit.Test;

import cohort.PercentileTransformer;

public class PercentileTransformerTest
{

    private final double[] COHORT_VALUES = { 0, 0, 5, 10, 10, 10, 10, 10, 10, 10 };
    private PercentileTransformer TRANSFORMER;

    @Before
    public void setUp()
    {
        TRANSFORMER = new PercentileTransformer(10);
        TRANSFORMER.fit(COHORT_VALUES);
    }

    @Test
    public void canDedupCohortValuesDuringFit()
    {
        double[] expectedDedupPercentiles = { 5, 20, 30, 70 };
        double[] expectedDedupRefValues = { 0, 4, 8.5, 10 };

        double[] actualDedupPercentiles = TRANSFORMER.getPercentilesDeduped();
        double[] actualDedupRefValues = TRANSFORMER.getRefValuesDeduped();

        assertArrayEquals(expectedDedupPercentiles, actualDedupPercentiles, 0.001);
        assertArrayEquals(expectedDedupRefValues, actualDedupRefValues, 0.001);
    }

    @Test
    public void canTransformInputValuesToPercentiles()
    {
        // Matches ref data exactly
        assertEquals(5,  TRANSFORMER.transform(0),   0.01);
        assertEquals(20, TRANSFORMER.transform(4),   0.01);
        assertEquals(30, TRANSFORMER.transform(8.5), 0.01);
        assertEquals(70, TRANSFORMER.transform(10),  0.01);

        // Requires interpolation
        assertEquals(12.50, TRANSFORMER.transform(2),  0.01);
        assertEquals(26.67, TRANSFORMER.transform(7),  0.01);

        // Out of bounds
        assertEquals(Double.POSITIVE_INFINITY, TRANSFORMER.transform( 1000),  0.01);
        assertEquals(Double.NEGATIVE_INFINITY, TRANSFORMER.transform(-1000),  0.01);
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
