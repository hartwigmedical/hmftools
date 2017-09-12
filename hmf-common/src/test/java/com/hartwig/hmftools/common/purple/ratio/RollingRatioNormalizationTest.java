package com.hartwig.hmftools.common.purple.ratio;

import static org.junit.Assert.assertEquals;

import java.util.List;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.purple.PurpleDatamodelTest;

import org.jetbrains.annotations.NotNull;
import org.junit.Test;

public class RollingRatioNormalizationTest {

    private static final String CHROMOSOME = "1";
    private static final double EPSILON = 1e-10;

    @Test
    public void testUsage() {
        final List<ReadRatio> input = Lists.newArrayList(
                create(1, 1.0),
                create(50, 1.5),
                create(120, 1.1),
                create(200, 1.2));

        final List<ReadRatio> output = new RollingRatioNormalization(100, input).get();
        assertEquals(input.size(), output.size());
        assertRatio(input.get(0), output.get(0), 1.25);
        assertRatio(input.get(1), output.get(1), 1.1);
        assertRatio(input.get(2), output.get(2), 1.2);
        assertRatio(input.get(3), output.get(3), 1.15);
    }

    private void assertRatio(@NotNull final ReadRatio input, @NotNull final ReadRatio output, double median) {
        assertEquals(input.position(), output.position());
        assertEquals(input.ratio() / median, output.ratio(), EPSILON);
    }

    private static ReadRatio create(long position, double ratio) {
        return PurpleDatamodelTest.createReadRatio(CHROMOSOME, position, ratio).build();
    }

}
