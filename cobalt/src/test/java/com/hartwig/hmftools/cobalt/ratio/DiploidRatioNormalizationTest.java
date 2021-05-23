package com.hartwig.hmftools.cobalt.ratio;

import static org.junit.Assert.assertEquals;

import java.util.List;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.cobalt.ImmutableReadRatio;
import com.hartwig.hmftools.common.cobalt.ReadRatio;

import org.jetbrains.annotations.NotNull;
import org.junit.Test;

public class DiploidRatioNormalizationTest
{

    private static final String CHROMOSOME = "1";
    private static final double EPSILON = 1e-10;

    @Test
    public void testCloseToZero()
    {
        final List<ReadRatio> input = Lists.newArrayList(create(0, 0), create(50, 0), create(100, 0.002), create(150, 0), create(200, 0));

        final List<ReadRatio> output = new DiploidRatioNormalization(1.0, 5, 5, input).get();
        assertEquals(input.size(), output.size());
        assertRatio(input.get(0), output.get(0), 1);
        assertRatio(input.get(1), output.get(1), 1);
        assertRatio(input.get(2), output.get(2), 1);
        assertRatio(input.get(3), output.get(3), 1);
        assertRatio(input.get(4), output.get(4), 1);
    }

    @Test
    public void testMaxWindowDistance()
    {
        final List<ReadRatio> input =
                Lists.newArrayList(create(0, 1.0), create(50, 1.5), create(100, -1), create(150, 1.1), create(200, 1.2));

        final List<ReadRatio> output = new DiploidRatioNormalization(1.0, 2, 1, input).get();
        assertEquals(input.size(), output.size());
        assertRatio(input.get(0), output.get(0), 1.25);
        assertRatio(input.get(1), output.get(1), 1.1);
        assertRatio(input.get(2), output.get(2), 1.0);
        assertRatio(input.get(3), output.get(3), 1.2);
        assertRatio(input.get(4), output.get(4), 1.15);
    }

    @Test
    public void testMinCoverage()
    {
        final List<ReadRatio> input =
                Lists.newArrayList(create(0, 1.0), create(50, 1.5), create(100, 2.0), create(150, -1), create(200, -1));

        final List<ReadRatio> output = new DiploidRatioNormalization(1.0, 1, 3, input).get();
        assertEquals(input.size(), output.size());
        assertRatio(input.get(0), output.get(0), 1.0);
        assertRatio(input.get(1), output.get(1), 1.5);
        assertRatio(input.get(2), output.get(2), 1.0);
        assertRatio(input.get(3), output.get(3), 1.0);
        assertRatio(input.get(4), output.get(4), 1.0);
    }

    private static void assertRatio(@NotNull final ReadRatio input, @NotNull final ReadRatio output, double median)
    {
        assertEquals(input.position(), output.position());
        assertEquals(input.ratio() / median, output.ratio(), EPSILON);
    }

    @NotNull
    private static ReadRatio create(long position, double ratio)
    {
        return createReadRatio(CHROMOSOME, position, ratio).build();
    }

    @NotNull
    private static ImmutableReadRatio.Builder createReadRatio(@NotNull final String chromosome, final long position, final double ratio)
    {
        return ImmutableReadRatio.builder().chromosome(chromosome).position(position).ratio(ratio);
    }
}
