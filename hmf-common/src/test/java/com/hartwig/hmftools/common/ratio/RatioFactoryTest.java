package com.hartwig.hmftools.common.ratio;

import com.hartwig.hmftools.common.exception.HartwigException;
import org.junit.Test;

import static org.junit.Assert.assertEquals;

public class RatioFactoryTest {

    private static final double EPSILON = 1.0e-10;

    @Test
    public void canConvertLine() throws HartwigException {
        final String line = "1\t1000\t0.3\t0.28\t2\t-1\t2\t-\t-1";

        final Ratio ratio = RatioFactory.fromRatioLine(line);
        assertEquals("1", ratio.chromosome());
        assertEquals(1001, ratio.position());
        assertEquals(0.3, ratio.ratio(), EPSILON);
        assertEquals(0.28, ratio.medianRatio(), EPSILON);
    }

}
