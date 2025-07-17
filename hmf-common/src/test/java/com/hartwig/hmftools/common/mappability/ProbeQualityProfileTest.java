package com.hartwig.hmftools.common.mappability;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertTrue;

import java.util.List;
import java.util.Map;
import java.util.OptionalDouble;

import com.google.common.io.Resources;
import com.hartwig.hmftools.common.region.ChrBaseRegion;

import org.junit.Test;

public class ProbeQualityProfileTest
{
    private final ProbeQualityProfile mProfile;

    public ProbeQualityProfileTest()
    {
        String testFile = Resources.getResource("mappability/test_probe_quality.tsv").getPath();
        mProfile = new ProbeQualityProfile(testFile);
    }

    @Test
    public void testConstructor()
    {
        Map<String, List<ProbeQualityWindow>> expectedWindows = Map.of(
                // General semi-realistic data.
                "1", List.of(
                        new ProbeQualityWindow(101, 141, 0.1f),
                        new ProbeQualityWindow(121, 161, 0.2f),
                        new ProbeQualityWindow(141, 181, 0.3f),
                        // Gap
                        new ProbeQualityWindow(181, 221, 0.4f),
                        new ProbeQualityWindow(201, 241, 0.5f)
                ),
                "2", List.of(
                        // No gaps
                        new ProbeQualityWindow(221, 261, 0.1f),
                        new ProbeQualityWindow(241, 281, 0.2f),
                        new ProbeQualityWindow(261, 301, 0.3f),
                        new ProbeQualityWindow(281, 321, 0.4f)
                ),
                "10", List.of(
                        new ProbeQualityWindow(761, 801, 0.1f),
                        new ProbeQualityWindow(781, 821, 0.2f),
                        // Gap
                        new ProbeQualityWindow(841, 881, 0.3f),
                        // Gap
                        new ProbeQualityWindow(901, 941, 0.4f),
                        new ProbeQualityWindow(921, 961, 0.5f),
                        // Gap
                        new ProbeQualityWindow(981, 1021, 0.6f)
                ),
                "11", List.of(
                        new ProbeQualityWindow(1, 41, 1.0f)
                ),
                // Data specifically for testing probe quality calculations.
                "LowQuality", List.of(
                        new ProbeQualityWindow(81, 121, 0.01f),
                        new ProbeQualityWindow(101, 141, 0.01f),
                        new ProbeQualityWindow(121, 161, 0.01f),
                        new ProbeQualityWindow(141, 181, 0.01f),
                        new ProbeQualityWindow(161, 201, 0.01f),
                        new ProbeQualityWindow(181, 221, 0.01f),
                        new ProbeQualityWindow(201, 241, 0.01f),
                        new ProbeQualityWindow(221, 261, 0.01f),
                        new ProbeQualityWindow(241, 281, 0.01f)
                ),
                "HighQuality", List.of(
                        new ProbeQualityWindow(81, 121, 1.0f),
                        new ProbeQualityWindow(101, 141, 1.0f),
                        new ProbeQualityWindow(121, 161, 1.0f),
                        new ProbeQualityWindow(141, 181, 1.0f),
                        new ProbeQualityWindow(161, 201, 1.0f),
                        new ProbeQualityWindow(181, 221, 1.0f),
                        new ProbeQualityWindow(201, 241, 1.0f),
                        new ProbeQualityWindow(221, 261, 1.0f),
                        new ProbeQualityWindow(241, 281, 1.0f)
                ),
                "MixedQuality", List.of(
                        new ProbeQualityWindow(81, 121, 0.1f),
                        new ProbeQualityWindow(101, 141, 0.2f),
                        new ProbeQualityWindow(121, 161, 0.4f),
                        new ProbeQualityWindow(141, 181, 0.3f),
                        new ProbeQualityWindow(161, 201, 0.6f),
                        new ProbeQualityWindow(181, 221, 0.5f),
                        new ProbeQualityWindow(201, 241, 0.8f),
                        new ProbeQualityWindow(221, 261, 0.7f),
                        new ProbeQualityWindow(241, 281, 0.9f)
                ),
                "QualityStep", List.of(
                        new ProbeQualityWindow(81, 121, 0.1f),
                        // 101-121 overlaps low and high quality.
                        new ProbeQualityWindow(101, 141, 1.0f),
                        new ProbeQualityWindow(121, 161, 1.0f),
                        new ProbeQualityWindow(141, 181, 1.0f),
                        new ProbeQualityWindow(161, 201, 1.0f),
                        new ProbeQualityWindow(181, 221, 1.0f),
                        new ProbeQualityWindow(201, 241, 1.0f),
                        new ProbeQualityWindow(221, 261, 1.0f),
                        // 241-261 overlaps low and high quality.
                        new ProbeQualityWindow(241, 281, 0.1f)
                )
        );
        assertEquals(expectedWindows, mProfile.mWindows);
    }

    @Test
    public void testComputeQualityScoreAllLow()
    {
        OptionalDouble actual = mProfile.computeQualityScore(new ChrBaseRegion("LowQuality", 101, 201));
        assertTrue(actual.isPresent());
        assertEquals(0.01, actual.getAsDouble(), 1e-6);
    }

    @Test
    public void testComputeQualityScoreAllHigh()
    {
        OptionalDouble actual = mProfile.computeQualityScore(new ChrBaseRegion("HighQuality", 101, 201));
        assertTrue(actual.isPresent());
        assertEquals(1.0, actual.getAsDouble(), 1e-6);
        // Quality score should never exceed 1.
        assertTrue(actual.getAsDouble() <= 1.0);
    }

    @Test
    public void testComputeQualityScoreBounds()
    {
        OptionalDouble actual = mProfile.computeQualityScore(new ChrBaseRegion("MixedQuality", 101, 201));
        assertTrue(actual.isPresent());
        // No need to duplicate the maths to find the exact score. It should be somewhere in between the scores of the windows.
        assertTrue(actual.getAsDouble() > 0.1 && actual.getAsDouble() < 0.8);
    }

    @Test
    public void testComputeQualityScoreEdge1()
    {
        OptionalDouble actual = mProfile.computeQualityScore(new ChrBaseRegion("QualityStep", 101, 201));
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
        OptionalDouble actual = mProfile.computeQualityScore(new ChrBaseRegion("QualityStep", 111, 201));
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
