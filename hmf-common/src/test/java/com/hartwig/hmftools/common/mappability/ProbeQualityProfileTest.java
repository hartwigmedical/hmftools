package com.hartwig.hmftools.common.mappability;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertTrue;

import java.util.List;
import java.util.Map;
import java.util.Optional;

import com.google.common.io.Resources;
import com.hartwig.hmftools.common.region.BaseRegion;
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
    public void test_constructor()
    {
        Map<String, List<ProbeQualityWindow>> expectedWindows = Map.of(
                "1", List.of(
                        new ProbeQualityWindow(101, 141, 0.1f),
                        new ProbeQualityWindow(121, 161, 0.2f),
                        new ProbeQualityWindow(141, 181, 0.3f),
                        new ProbeQualityWindow(181, 221, 0.4f),
                        new ProbeQualityWindow(201, 241, 0.5f)
                ),
                "2", List.of(
                        new ProbeQualityWindow(221, 261, 0.1f),
                        new ProbeQualityWindow(241, 281, 0.2f),
                        new ProbeQualityWindow(261, 301, 0.3f),
                        new ProbeQualityWindow(281, 321, 0.4f)
                ),
                "10", List.of(
                        new ProbeQualityWindow(761, 801, 0.1f),
                        new ProbeQualityWindow(781, 821, 0.2f),
                        new ProbeQualityWindow(841, 881, 0.3f),
                        new ProbeQualityWindow(901, 941, 0.4f),
                        new ProbeQualityWindow(921, 961, 0.5f),
                        new ProbeQualityWindow(981, 1021, 0.6f)
                ),
                "11", List.of(
                        new ProbeQualityWindow(1, 41, 1.0f)
                )
        );
        Map<String, List<BaseRegion>> expectedCoverage = Map.of(
                "1", List.of(
                        new BaseRegion(101, 241)
                ),
                "2", List.of(
                        new BaseRegion(221, 321)
                ),
                "10", List.of(
                        new BaseRegion(761, 821),
                        new BaseRegion(841, 881),
                        new BaseRegion(901, 961),
                        new BaseRegion(981, 1021)
                ),
                "11", List.of(
                        new BaseRegion(1, 41)
                )
        );
        assertEquals(expectedWindows, mProfile.mWindows);
        assertEquals(expectedCoverage, mProfile.mCoverage);
    }

    @Test
    public void test_computeQualityScore_covered()
    {
        Optional<Double> actual;

        // Probe contained within 1 window.
        actual = mProfile.computeQualityScore(new ChrBaseRegion("1", 101, 111));
        assertTrue(actual.isPresent());
        assertEquals(0.1, actual.get(), 1e-6);

        // Probe overlapping 2 windows.
        actual = mProfile.computeQualityScore(new ChrBaseRegion("1", 101, 131));
        assertTrue(actual.isPresent());
        // No need to duplicate the maths to find the exact score. It should be somewhere in between the scores of the 2 windows.
        assertTrue(actual.get() > 0.1 && actual.get() < 0.2);

        // Probe overlapping multiple windows.
        actual = mProfile.computeQualityScore(new ChrBaseRegion("2", 233, 282));
        assertTrue(actual.isPresent());
        assertTrue(actual.get() > 0.1 && actual.get() < 0.3);
    }

    @Test
    public void test_computeQualityScore_notCovered()
    {
        // Nonexistent chromosome.
        assertEquals(Optional.empty(), mProfile.computeQualityScore(new ChrBaseRegion("50", 30, 50)));
        // Region not covered at all.
        assertEquals(Optional.empty(), mProfile.computeQualityScore(new ChrBaseRegion("1", 1, 15)));
        assertEquals(Optional.empty(), mProfile.computeQualityScore(new ChrBaseRegion("1", 10000, 20000)));
        // Region partially covered.
        assertEquals(Optional.empty(), mProfile.computeQualityScore(new ChrBaseRegion("1", 50, 160)));
        assertEquals(Optional.empty(), mProfile.computeQualityScore(new ChrBaseRegion("1", 145, 1000)));
        assertEquals(Optional.empty(), mProfile.computeQualityScore(new ChrBaseRegion("1", 50, 1000)));
    }
}
