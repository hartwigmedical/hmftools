package com.hartwig.hmftools.compar.common.field;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertFalse;
import static org.junit.Assert.assertSame;
import static org.junit.Assert.assertTrue;

import java.util.HashSet;
import java.util.List;
import java.util.Set;

import com.hartwig.hmftools.common.gene.TranscriptRegionType;
import com.hartwig.hmftools.compar.linx.BreakendData;
import com.hartwig.hmftools.compar.linx.TestBreakendDataBuilder;

import org.junit.Test;

public class BreakendsFieldTest
{
    private static final String FIELD_NAME = "Breakends";

    private static BreakendsField field(final boolean isCompared)
    {
        return new BreakendsField(FIELD_NAME, i -> ((TestFieldItem<List<BreakendData>>) i).Value, isCompared);
    }

    @Test
    public void nameAndIsComparedReflectConstructorArgs()
    {
        assertEquals(FIELD_NAME, field(true).name());
        assertTrue(field(true).isCompared());
        assertFalse(field(false).isCompared());
    }

    @Test
    public void displayValueJoinsFullStringsWithDelimiter()
    {
        BreakendData breakend1 = TestBreakendDataBuilder.BUILDER.create();
        BreakendData breakend2 = TestBreakendDataBuilder.BUILDER.createWithAlternateDefaults();

        String expected = breakend1.fullStr(true) + Field.DISPLAY_VALUE_DELIMITER + breakend2.fullStr(true);
        assertEquals(expected, field(true).displayValue(new TestFieldItem<>(List.of(breakend1, breakend2))));
    }

    @Test
    public void hasDiffIsFalseWhenBreakendsMatchAndHaveEqualAttributes()
    {
        BreakendData oldBreakend = TestBreakendDataBuilder.BUILDER.create();
        BreakendData newBreakend = TestBreakendDataBuilder.BUILDER.create();

        assertFalse(field(true).hasDiff(new TestFieldItem<>(List.of(oldBreakend)), new TestFieldItem<>(List.of(newBreakend))));
    }

    @Test
    public void hasDiffIsTrueWhenABreakendIsUnmatched()
    {
        BreakendData oldBreakend = TestBreakendDataBuilder.BUILDER.create();
        BreakendData newBreakend = TestBreakendDataBuilder.BUILDER.create(b ->
        {
            b.position = 1000000;
            b.comparisonPosition = 1000000;
        });

        assertTrue(field(true).hasDiff(new TestFieldItem<>(List.of(oldBreakend)), new TestFieldItem<>(List.of(newBreakend))));
    }

    @Test
    public void hasDiffIsTrueWhenMatchedBreakendsHaveDifferentReportedStatus()
    {
        BreakendData oldBreakend = TestBreakendDataBuilder.BUILDER.create();
        BreakendData newBreakend = TestBreakendDataBuilder.BUILDER.create(b -> b.reportedDisruption = false);

        assertTrue(field(true).hasDiff(new TestFieldItem<>(List.of(oldBreakend)), new TestFieldItem<>(List.of(newBreakend))));
    }

    @Test
    public void determineDiffsIsEmptyWhenNoDiff()
    {
        BreakendData oldBreakend = TestBreakendDataBuilder.BUILDER.create();
        BreakendData newBreakend = TestBreakendDataBuilder.BUILDER.create();

        List<String> diffs =
                field(true).determineDiffs(new TestFieldItem<>(List.of(oldBreakend)), new TestFieldItem<>(List.of(newBreakend)));
        assertTrue(diffs.isEmpty());
    }

    @Test
    public void determineDiffsReportsUnmatchedBreakendsOnBothSides()
    {
        BreakendData oldBreakend = TestBreakendDataBuilder.BUILDER.create();
        BreakendData newBreakend = TestBreakendDataBuilder.BUILDER.create(b ->
        {
            b.position = 1000000;
            b.comparisonPosition = 1000000;
        });

        List<String> diffs =
                field(true).determineDiffs(new TestFieldItem<>(List.of(oldBreakend)), new TestFieldItem<>(List.of(newBreakend)));

        Set<String> expected = Set.of(
                String.format("%s(unmatched SV %s/)", FIELD_NAME, oldBreakend.svInfoStr()),
                String.format("%s(unmatched SV /%s)", FIELD_NAME, newBreakend.svInfoStr()));
        assertEquals(expected, new HashSet<>(diffs));
    }

    @Test
    public void determineDiffsReportsTranscriptDetailMismatchForMatchedBreakends()
    {
        BreakendData oldBreakend = TestBreakendDataBuilder.BUILDER.create();
        BreakendData newBreakend = TestBreakendDataBuilder.BUILDER.create(b -> b.regionType = TranscriptRegionType.EXONIC);

        List<String> diffs =
                field(true).determineDiffs(new TestFieldItem<>(List.of(oldBreakend)), new TestFieldItem<>(List.of(newBreakend)));

        String expected = String.format("%s(%s transcript %s/%s)",
                FIELD_NAME, oldBreakend.svInfoStr(), oldBreakend.transcriptStr(), newBreakend.transcriptStr());
        assertEquals(List.of(expected), diffs);
    }

    @Test
    public void determineDiffsReportsReportedStatusMismatchForMatchedBreakends()
    {
        BreakendData oldBreakend = TestBreakendDataBuilder.BUILDER.create();
        BreakendData newBreakend = TestBreakendDataBuilder.BUILDER.create(b -> b.reportedDisruption = false);

        List<String> diffs =
                field(true).determineDiffs(new TestFieldItem<>(List.of(oldBreakend)), new TestFieldItem<>(List.of(newBreakend)));

        String expected = String.format("%s(%s reported %s/%s)",
                FIELD_NAME, oldBreakend.svInfoStr(), oldBreakend.Breakend.reportedStatus(), newBreakend.Breakend.reportedStatus());
        assertEquals(List.of(expected), diffs);
    }

    @Test
    public void withComparedReturnsNewFieldWithUpdatedIsCompared()
    {
        BreakendsField field = field(true);
        Field updated = field.withCompared(false);
        assertFalse(updated.isCompared());
    }

    @Test
    public void withThresholdsAreUnsupportedAndReturnSameInstance()
    {
        BreakendsField field = field(true);
        assertSame(field, field.withAbsoluteThreshold(5.0));
        assertSame(field, field.withPercentThreshold(0.2));
    }
}
