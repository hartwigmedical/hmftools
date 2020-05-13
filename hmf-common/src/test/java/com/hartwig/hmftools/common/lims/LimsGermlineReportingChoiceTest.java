package com.hartwig.hmftools.common.lims;

import static org.junit.Assert.assertEquals;

import org.junit.Test;

public class LimsGermlineReportingChoiceTest {

    @Test
    public void canExtractGermlineLevel() {
        assertEquals(LimsGermlineReportingLevel.REPORT_WITH_NOTIFICATION,
                LimsGermlineReportingLevel.fromLimsGermlineReportingLevelString("1: Behandelbare toevalsbevindingen", "WIDE02991111T", true));
        assertEquals(LimsGermlineReportingLevel.REPORT_WITH_NOTIFICATION,
                LimsGermlineReportingLevel.fromLimsGermlineReportingLevelString("2: Alle toevalsbevindingen", "WIDE02991111T", true));
        assertEquals(LimsGermlineReportingLevel.REPORT_WITHOUT_NOTIFICATION,
                LimsGermlineReportingLevel.fromLimsGermlineReportingLevelString(
                        "3: Geen toevalsbevindingen; familie mag deze wel opvragen",
                        "WIDE02991111T", true));
        assertEquals(LimsGermlineReportingLevel.REPORT_WITHOUT_NOTIFICATION,
                LimsGermlineReportingLevel.fromLimsGermlineReportingLevelString("3: Geen toevalsbevindingen", "WIDE02991111T", true));
        assertEquals(LimsGermlineReportingLevel.REPORT_WITHOUT_NOTIFICATION,
                LimsGermlineReportingLevel.fromLimsGermlineReportingLevelString(
                        "4: Geen toevalsbevindingen; familie mag deze niet opvragen",
                        "WIDE02991111T", true));
        assertEquals(LimsGermlineReportingLevel.NO_REPORTING,
                LimsGermlineReportingLevel.fromLimsGermlineReportingLevelString("", "CPCT02991111T", false));
        assertEquals(LimsGermlineReportingLevel.NO_REPORTING,
                LimsGermlineReportingLevel.fromLimsGermlineReportingLevelString("", "DRUP02991111T", false));
        assertEquals(LimsGermlineReportingLevel.NO_REPORTING,
                LimsGermlineReportingLevel.fromLimsGermlineReportingLevelString("", "COLO02991111T", false));
        assertEquals(LimsGermlineReportingLevel.REPORT_WITHOUT_NOTIFICATION,
                LimsGermlineReportingLevel.fromLimsGermlineReportingLevelString("", "CORE02991111T", true));

        assertEquals(LimsGermlineReportingLevel.NO_REPORTING,
                LimsGermlineReportingLevel.fromLimsGermlineReportingLevelString("", "CORE02991111T", false));
        assertEquals(LimsGermlineReportingLevel.NO_REPORTING,
                LimsGermlineReportingLevel.fromLimsGermlineReportingLevelString("3: Geen toevalsbevindingen", "WIDE02991111T", false));
    }

    @Test(expected = IllegalStateException.class)
    public void hasUnknownGermlineChoice() {
        LimsGermlineReportingLevel.fromLimsGermlineReportingLevelString("ALL", "WIDE02991111T", true);
    }
}