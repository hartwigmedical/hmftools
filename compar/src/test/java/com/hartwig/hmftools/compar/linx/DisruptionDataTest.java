package com.hartwig.hmftools.compar.linx;

import static com.hartwig.hmftools.compar.linx.DisruptionData.FLD_UNMATCHED_SV;

import static org.junit.Assert.assertNull;
import static org.junit.Assert.assertTrue;

import static junit.framework.TestCase.assertEquals;

import java.util.Collections;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.function.Consumer;

import com.hartwig.hmftools.compar.ComparConfig;
import com.hartwig.hmftools.compar.ComparableItemTest;
import com.hartwig.hmftools.compar.common.DiffThresholds;
import com.hartwig.hmftools.compar.common.MatchLevel;
import com.hartwig.hmftools.compar.common.Mismatch;
import com.hartwig.hmftools.compar.common.MismatchType;

import org.junit.Before;
import org.junit.Test;

public class DisruptionDataTest extends ComparableItemTest<DisruptionData, DisruptionComparer, TestDisruptionDataBuilder>
{
    @Before
    public void setUp()
    {
        comparer = new DisruptionComparer(new ComparConfig());
        builder = TestDisruptionDataBuilder.BUILDER;

        DisruptionData alternateValueSource = builder.createWithAlternateDefaults();

        // Does not include every field because field comparisons within breakends don't work well in generic tests
        fieldToAlternateValueInitializer = Map.of(FLD_UNMATCHED_SV, b -> b.breakends = alternateValueSource.Breakends);

        nameToAlternateIndexInitializer = Map.of("GeneName", b -> b.geneName = alternateValueSource.GeneName);
        reportabilityFieldToFalseReportabilityInitializer = Collections.emptyMap();
    }

    @Test
    public void fullyMatchesSelfWithLiftover()
    {
        BreakendData refBreakend = TestBreakendDataBuilder.BUILDER.create(b -> {
            b.comparisonChromosome = "21";
            b.comparisonPosition = 40000000;
        });
        BreakendData newBreakend = TestBreakendDataBuilder.BUILDER.create(b -> {
            b.chromosome = "21";
            b.position = 40000000;
        });
        DisruptionData refVictim = builder.create(b -> b.breakends = List.of(refBreakend));
        DisruptionData newVictim = builder.create(b -> b.breakends = List.of(newBreakend));

        DiffThresholds diffThresholds = createDefaultThresholds();

        assertTrue(refVictim.matches(newVictim));
        assertNull(refVictim.findMismatch(newVictim, MatchLevel.DETAILED, diffThresholds, false));
    }

    @Test
    public void unmatchedSvsAreRecognized()
    {
        BreakendData sharedBreakend1 = TestBreakendDataBuilder.BUILDER.create();
        BreakendData sharedBreakend2 = TestBreakendDataBuilder.BUILDER.create(b -> {
            b.position = 41100000;
            b.comparisonPosition = 41100000;
        });
        BreakendData refUniqueBreakend = TestBreakendDataBuilder.BUILDER.create(b -> {
            b.position = 42000000;
            b.comparisonPosition = 42000000;
        });
        BreakendData newUniqueBreakend = TestBreakendDataBuilder.BUILDER.create(b -> {
            b.position = 40000000;
            b.comparisonPosition = 40000000;
        });
        DisruptionData refVictim = builder.create(b -> b.breakends = List.of(sharedBreakend1, sharedBreakend2, refUniqueBreakend));
        DisruptionData newVictim = builder.create(b -> b.breakends = List.of(newUniqueBreakend, sharedBreakend1, sharedBreakend2));

        DiffThresholds diffThresholds = createDefaultThresholds();

        assertTrue(refVictim.matches(newVictim));

        Mismatch mismatch = refVictim.findMismatch(newVictim, MatchLevel.DETAILED, diffThresholds, false);

        assertEquals(MismatchType.VALUE, mismatch.MismatchType());
        assertEquals(refVictim, mismatch.RefItem());
        assertEquals(newVictim, mismatch.NewItem());
        Set<String> expectedDiffs = Set.of("unmatchedSv(:BND chr21:42000000:0/)", "unmatchedSv(/:BND chr21:40000000:0)");
        assertEquals(expectedDiffs, new HashSet<>(mismatch.DiffValues()));
    }

    @Test
    public void regionTypeMismatchInBreakendIsRecognized()
    {
        BreakendData alternateValueSource = TestBreakendDataBuilder.BUILDER.createWithAlternateDefaults();
        Consumer<TestBreakendDataBuilder> altInitializer = b -> b.regionType = alternateValueSource.Breakend.regionType();
        String field = "RegionType";
        String expectedDiff =
                "breakend(:BND chr21:41500000:0 transcript ENST00000332149:CODING:INTRONIC:-1/ENST00000332149:CODING:EXONIC:-1)";

        assertSingleFieldDifferenceInBreakendRecognized(field, altInitializer, expectedDiff);
    }

    @Test
    public void codingTypeMismatchInBreakendIsRecognized()
    {
        BreakendData alternateValueSource = TestBreakendDataBuilder.BUILDER.createWithAlternateDefaults();
        Consumer<TestBreakendDataBuilder> altInitializer = b -> b.codingType = alternateValueSource.Breakend.codingType();
        String field = "CodingType";
        String expectedDiff =
                "breakend(:BND chr21:41500000:0 transcript ENST00000332149:CODING:INTRONIC:-1/ENST00000332149:NON_CODING:INTRONIC:-1)";

        assertSingleFieldDifferenceInBreakendRecognized(field, altInitializer, expectedDiff);
    }

    @Test
    public void nextSpliceExonRankMismatchInBreakendIsRecognized()
    {
        BreakendData alternateValueSource = TestBreakendDataBuilder.BUILDER.createWithAlternateDefaults();
        Consumer<TestBreakendDataBuilder> altInitializer = b -> b.nextSpliceExonRank = alternateValueSource.Breakend.nextSpliceExonRank();
        String field = "NextSpliceExonRank";
        String expectedDiff =
                "breakend(:BND chr21:41500000:0 transcript ENST00000332149:CODING:INTRONIC:-1/ENST00000332149:CODING:INTRONIC:-1)";

        assertSingleFieldDifferenceInBreakendRecognized(field, altInitializer, expectedDiff);
    }

    @Test
    public void reportabilityDifferenceIsRecognizedInDetailedMode()
    {
        DiffThresholds diffThresholds = createDefaultThresholds();
        BreakendData defaultBreakend = TestBreakendDataBuilder.BUILDER.create();
        BreakendData nonDefaultBreakend = TestBreakendDataBuilder.BUILDER.create(b -> b.reportedDisruption = false);

        DisruptionData refVictim = builder.create(b -> b.breakends = List.of(defaultBreakend));
        DisruptionData newVictim = builder.create(b -> b.breakends = List.of(nonDefaultBreakend));

        assertTrue(refVictim.matches(newVictim));
        Mismatch mismatch = refVictim.findMismatch(newVictim, MatchLevel.DETAILED, diffThresholds, false);

        assertEquals(MismatchType.VALUE, mismatch.MismatchType());
        assertEquals(refVictim, mismatch.RefItem());
        assertEquals(newVictim, mismatch.NewItem());
        assertEquals(List.of("breakend(:BND chr21:41500000:0 reported REPORTED/NONE)"), mismatch.DiffValues());
    }

    @Test
    public void reportabilityDifferenceIsRecognizedInReportableMode()
    {
        DiffThresholds diffThresholds = createDefaultThresholds();
        BreakendData defaultBreakend = TestBreakendDataBuilder.BUILDER.create();
        BreakendData nonDefaultBreakend = TestBreakendDataBuilder.BUILDER.create(b -> b.reportedDisruption = false);

        DisruptionData reportableVictim = builder.create(b -> b.breakends = List.of(defaultBreakend));
        DisruptionData nonReportableVictim = builder.create(b -> b.breakends = List.of(nonDefaultBreakend));

        assertTrue(reportableVictim.matches(nonReportableVictim));
        Mismatch mismatch = reportableVictim.findMismatch(nonReportableVictim, MatchLevel.REPORTABLE, diffThresholds, false);

        assertEquals(MismatchType.REF_ONLY, mismatch.MismatchType());
        assertEquals(reportableVictim, mismatch.RefItem());
        assertEquals(nonReportableVictim, mismatch.NewItem());
        assertEquals(List.of("breakend(:BND chr21:41500000:0 reported REPORTED/NONE)"), mismatch.DiffValues());

        Mismatch oppositeMismatch = nonReportableVictim.findMismatch(reportableVictim, MatchLevel.REPORTABLE, diffThresholds, false);

        assertEquals(MismatchType.NEW_ONLY, oppositeMismatch.MismatchType());
        assertEquals(nonReportableVictim, oppositeMismatch.RefItem());
        assertEquals(reportableVictim, oppositeMismatch.NewItem());
        assertEquals(List.of("breakend(:BND chr21:41500000:0 reported NONE/REPORTED)"), oppositeMismatch.DiffValues());
    }

    private void assertSingleFieldDifferenceInBreakendRecognized(final String field, final Consumer<TestBreakendDataBuilder> altInitializer,
            final String expectedDiff)
    {
        DiffThresholds diffThresholds = createDefaultThresholds();

        BreakendData defaultBreakend = TestBreakendDataBuilder.BUILDER.create();
        BreakendData nonDefaultBreakend = TestBreakendDataBuilder.BUILDER.create(altInitializer);

        DisruptionData refVictim = builder.create(b -> b.breakends = List.of(defaultBreakend));
        DisruptionData newVictim = builder.create(b -> b.breakends = List.of(nonDefaultBreakend));

        assertTrue("Test difference in " + field, refVictim.matches(newVictim));

        Mismatch mismatch = refVictim.findMismatch(newVictim, MatchLevel.DETAILED, diffThresholds, false);

        assertEquals("Test difference in " + field, MismatchType.VALUE, mismatch.MismatchType());
        assertEquals("Test difference in " + field, refVictim, mismatch.RefItem());
        assertEquals("Test difference in " + field, newVictim, mismatch.NewItem());
        assertEquals("Test difference in " + field, List.of(expectedDiff), mismatch.DiffValues());
    }
}
