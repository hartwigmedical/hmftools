package com.hartwig.hmftools.common.mappability;

import static org.junit.Assert.assertArrayEquals;
import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertTrue;

import java.util.OptionalDouble;

import com.google.common.io.Resources;
import com.hartwig.hmftools.common.region.BaseRegion;
import com.hartwig.hmftools.common.region.ChrBaseRegion;

import org.junit.Test;

public class ProbeQualityProfileTest
{
    private static final int BASE_WINDOW_LENGTH = 40;
    private static final float EPSILON = 1e-9f;

    private final ProbeQualityProfile mProfile;

    public ProbeQualityProfileTest()
    {
        String testFile = Resources.getResource("mappability/test_probe_quality.tsv").getPath();
        mProfile = ProbeQualityProfile.loadFromResourceFile(testFile);
    }

    private static void assertWindowArrayEqual(final ProbeQualityProfile.WindowArray actual, final int[] expectedEndPositions,
            final float[] expectedQualityScores)
    {
        assertEquals(BASE_WINDOW_LENGTH, actual.mBaseWindowLength);
        assertEquals(expectedEndPositions.length, expectedQualityScores.length);
        assertEquals(expectedEndPositions.length, actual.size());
        assertArrayEquals(expectedEndPositions, actual.mEndPositions);
        assertArrayEquals(expectedQualityScores, actual.mQualityScores, EPSILON);
        for(int i = 0; i < expectedEndPositions.length; ++i)
        {
            assertEquals(actual.getRegion(i), new BaseRegion(expectedEndPositions[i] - BASE_WINDOW_LENGTH + 1, expectedEndPositions[i]));
            assertEquals(actual.getQualityScore(i), expectedQualityScores[i], EPSILON);
        }
    }

    @Test
    public void testConstructor()
    {
        // General semi-realistic data.
        assertWindowArrayEqual(mProfile.mWindows.get("1"),
                // Gap
                new int[] { 140, 160, 180, 220, 240 },
                new float[] { 0.1f, 0.2f, 0.3f, 0.4f, 0.5f });
        assertWindowArrayEqual(mProfile.mWindows.get("2"),
                // No gaps
                new int[] { 260, 280, 300, 320 },
                new float[] { 0.1f, 0.2f, 0.3f, 0.4f });
        assertWindowArrayEqual(mProfile.mWindows.get("10"),
                // Multiple gaps
                new int[] { 800, 820, 880, 940, 960, 1020 },
                new float[] { 0.1f, 0.2f, 0.3f, 0.4f, 0.5f, 0.6f });
        assertWindowArrayEqual(mProfile.mWindows.get("11"),
                // Single window
                new int[] { 40 },
                new float[] { 1.0f });
        // Data specifically for testing probe quality calculations.
        assertWindowArrayEqual(mProfile.mWindows.get("LowQuality"),
                new int[] { 120, 140, 160, 180, 200, 220, 240, 260, 280 },
                new float[] { 0.01f, 0.01f, 0.01f, 0.01f, 0.01f, 0.01f, 0.01f, 0.01f, 0.01f });
        assertWindowArrayEqual(mProfile.mWindows.get("HighQuality"),
                new int[] { 120, 140, 160, 180, 200, 220, 240, 260, 280 },
                new float[] { 1.0f, 1.0f, 1.0f, 1.0f, 1.0f, 1.0f, 1.0f, 1.0f, 1.0f });
        assertWindowArrayEqual(mProfile.mWindows.get("MixedQuality"),
                new int[] { 120, 140, 160, 180, 200, 220, 240, 260, 280 },
                new float[] { 0.1f, 0.2f, 0.4f, 0.3f, 0.6f, 0.5f, 0.8f, 0.7f, 0.9f });
        assertWindowArrayEqual(mProfile.mWindows.get("QualityStep"),
                // 101-120 overlaps low and high quality.
                // 241-260 overlaps low and high quality.
                new int[] { 120, 140, 160, 180, 200, 220, 240, 260, 280 },
                new float[] { 0.1f, 1.0f, 1.0f, 1.0f, 1.0f, 1.0f, 1.0f, 1.0f, 0.1f });
    }

    @Test
    public void testComputeQualityScoreAllLow()
    {
        OptionalDouble actual = mProfile.computeQualityScore(new ChrBaseRegion("LowQuality", 101, 200));
        assertTrue(actual.isPresent());
        assertEquals(0.01, actual.getAsDouble(), 1e-6);
    }

    @Test
    public void testComputeQualityScoreAllHigh()
    {
        OptionalDouble actual = mProfile.computeQualityScore(new ChrBaseRegion("HighQuality", 101, 200));
        assertTrue(actual.isPresent());
        assertEquals(1.0, actual.getAsDouble(), 1e-6);
        // Quality score should never exceed 1.
        assertTrue(actual.getAsDouble() <= 1.0);
    }

    @Test
    public void testComputeQualityScoreBounds()
    {
        OptionalDouble actual = mProfile.computeQualityScore(new ChrBaseRegion("MixedQuality", 101, 200));
        assertTrue(actual.isPresent());
        // No need to duplicate the maths to find the exact score. It should be somewhere in between the scores of the windows.
        assertTrue(actual.getAsDouble() > 0.1 && actual.getAsDouble() < 0.8);
    }

    @Test
    public void testComputeQualityScoreEdge1()
    {
        OptionalDouble actual = mProfile.computeQualityScore(new ChrBaseRegion("QualityStep", 101, 200));
        assertTrue(actual.isPresent());
        // Interpolating between 0.1 and 1.0 quality score, overlaps of 20b and 20b.
        // Risk1 ~= 1/0.1 - 1 = 9
        // Risk2 ~= 1/1.0 - 1 = 0
        // Risk ~= 9 * 20/40 + 0 * 20/40 = 4.5
        // Quality = 1 / (4.5 + 1) ~= 0.182
        assertEquals(0.182, actual.getAsDouble(), 0.1);
    }

    @Test
    public void testComputeQualityScoreEdge2()
    {
        OptionalDouble actual = mProfile.computeQualityScore(new ChrBaseRegion("QualityStep", 111, 200));
        assertTrue(actual.isPresent());
        // Interpolating between 0.1 and 1.0 quality score, with overlaps of 10 and 30
        // Risk1 ~= 1/0.1 - 1 = 9
        // Risk2 ~= 1/1.0 - 1 = 0
        // Risk ~= 9 * 10/40 + 0 * 30/40 = 2.25
        // Quality = 1 / (2.25 + 1) ~= 0.307
        assertEquals(0.307, actual.getAsDouble(), 0.1);
    }

    @Test
    public void testComputeQualityScoreNotCovered()
    {
        // Nonexistent chromosome.
        assertEquals(OptionalDouble.empty(), mProfile.computeQualityScore(new ChrBaseRegion("50", 30, 150)));
        // Region not covered at all.
        assertEquals(OptionalDouble.empty(), mProfile.computeQualityScore(new ChrBaseRegion("1", 1, 120)));
        assertEquals(OptionalDouble.empty(), mProfile.computeQualityScore(new ChrBaseRegion("1", 10000, 20000)));
        // Region partially covered.
        assertEquals(OptionalDouble.empty(), mProfile.computeQualityScore(new ChrBaseRegion("1", 20, 160)));
        assertEquals(OptionalDouble.empty(), mProfile.computeQualityScore(new ChrBaseRegion("1", 145, 1000)));
        assertEquals(OptionalDouble.empty(), mProfile.computeQualityScore(new ChrBaseRegion("1", 50, 1000)));
    }
}
