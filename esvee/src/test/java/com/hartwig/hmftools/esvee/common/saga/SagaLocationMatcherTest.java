package com.hartwig.hmftools.esvee.common.saga;

import static java.lang.String.format;

import static com.hartwig.hmftools.common.genome.region.Orientation.FORWARD;
import static com.hartwig.hmftools.common.genome.region.Orientation.REVERSE;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertNotNull;
import static org.junit.Assert.assertNull;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

import com.hartwig.hmftools.common.genome.region.Orientation;
import com.hartwig.hmftools.common.region.BasePosition;

import org.junit.Test;

public class SagaLocationMatcherTest
{
    private static final int MAX_DISTANCE = 50;
    private static final String CHROMOSOME = "chr1";

    private final SagaLocationMatcher mMatcher;
    private final Map<String, List<SagaIndexedBreakend>> mBreakendMap;
    private int mVariantIdCounter; // Auto-generate unique variant IDs

    public SagaLocationMatcherTest()
    {
        // Breakends used by the matcher, to be populated later in the test case.
        mBreakendMap = new HashMap<>();
        mVariantIdCounter = 0;
        mMatcher = new SagaLocationMatcher(
                new SagaLocationMatcher.Config(MAX_DISTANCE),
                mBreakendMap
        );
    }

    private void addBreakend(final String chromosome, int position, final Orientation orientation)
    {
        String variantId = format("var%d", ++mVariantIdCounter);

        // Create a minimal valid SagaVariant that satisfies all invariants.
        SagaVariant variant = new SagaVariant(
                variantId,
                new SagaBreakend(new BasePosition(chromosome, position), FORWARD),
                // The second breakend is not needed and not relevant, so put it very far away.
                new SagaBreakend(new BasePosition(chromosome, position + 1_000_000), REVERSE),
                ""
        );

        SagaIndexedBreakend breakend = new SagaIndexedBreakend(
                new SagaBreakend(new BasePosition(chromosome, position), orientation),
                variant
        );
        addBreakend(breakend);
    }

    private void addBreakend(final SagaIndexedBreakend breakend)
    {
        List<SagaIndexedBreakend> chrBreakends = mBreakendMap.computeIfAbsent(breakend.chromosome(), k -> new ArrayList<>());
        chrBreakends.add(breakend);
        chrBreakends.sort(null); // Sort by position (SagaIndexedBreakend implements Comparable)
    }

    @Test
    public void testNoSagaBreakendOnTargetChromosome()
    {
        addBreakend("chr2", 100, FORWARD);

        SagaLocationMatch result = mMatcher.match(CHROMOSOME, 100, FORWARD);

        assertNull("Should return null when no breakends exist on chromosome", result);
    }

    @Test
    public void testNoSagaBreakendWithinMaxDistance()
    {
        addBreakend(CHROMOSOME, 100, FORWARD);

        SagaLocationMatch result = mMatcher.match(CHROMOSOME, 151, FORWARD);

        assertNull("Should return null when nearest breakend exceeds max distance", result);
    }

    @Test
    public void testExactly1SagaBreakendWithinMaxDistance()
    {
        addBreakend(CHROMOSOME, 100, FORWARD);

        SagaLocationMatch result = mMatcher.match(CHROMOSOME, 120, FORWARD);

        assertNotNull("Should return a match for the single breakend within distance", result);
        assertEquals("var1", result.variantId());
        assertEquals(20, result.distance());
    }

    @Test
    public void testMultipleSagaBreakendsDifferentDistances()
    {
        addBreakend(CHROMOSOME, 80, FORWARD);
        addBreakend(CHROMOSOME, 110, FORWARD);
        addBreakend(CHROMOSOME, 145, FORWARD);

        SagaLocationMatch result = mMatcher.match(CHROMOSOME, 100, FORWARD);

        assertNotNull("Should return the closest breakend", result);
        assertEquals("var2", result.variantId());
        assertEquals(10, result.distance());
    }

    @Test
    public void testMultipleSagaBreakendsSameDistance()
    {
        addBreakend(CHROMOSOME, 95, FORWARD);
        addBreakend(CHROMOSOME, 105, REVERSE);

        SagaLocationMatch result = mMatcher.match(CHROMOSOME, 100, FORWARD);

        assertNotNull("Should return a match for the tie-break", result);
        assertEquals(5, result.distance());
        assertEquals("var1", result.variantId());
    }

    @Test
    public void testMultipleSagaBreakendsSameDistanceThenVariantId()
    {
        addBreakend(CHROMOSOME, 95, FORWARD);
        addBreakend(CHROMOSOME, 105, FORWARD);

        SagaLocationMatch result = mMatcher.match(CHROMOSOME, 100, FORWARD);

        assertNotNull("Should return a match", result);
        assertEquals(5, result.distance());
        assertEquals("var1", result.variantId());
    }

    @Test
    public void testSagaBreakendAtExactPosition()
    {
        addBreakend(CHROMOSOME, 100, FORWARD);

        SagaLocationMatch result = mMatcher.match(CHROMOSOME, 100, FORWARD);

        assertNotNull("Should return a match at exact position", result);
        assertEquals("var1", result.variantId());
        assertEquals(0, result.distance());
    }

    @Test
    public void testMultipleBreakendAtExactPosition()
    {
        addBreakend(CHROMOSOME, 100, FORWARD);
        addBreakend(CHROMOSOME, 100, FORWARD);

        SagaLocationMatch result = mMatcher.match(CHROMOSOME, 100, FORWARD);

        assertNotNull("Should return a match", result);
        assertEquals(0, result.distance());
        assertEquals("var1", result.variantId());
    }

    @Test
    public void testBoundaryAtMaxDistance()
    {
        addBreakend(CHROMOSOME, 100, FORWARD);

        SagaLocationMatch result = mMatcher.match(CHROMOSOME, 150, FORWARD);

        assertNotNull("Should return a match at exactly max distance", result);
        assertEquals("var1", result.variantId());
        assertEquals(50, result.distance());
    }

    @Test
    public void testJustBeyondMaxDistance()
    {
        addBreakend(CHROMOSOME, 100, FORWARD);

        SagaLocationMatch result = mMatcher.match(CHROMOSOME, 151, FORWARD);

        assertNull("Should return null when distance exceeds max distance by 1", result);
    }
}