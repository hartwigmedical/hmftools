package com.hartwig.hmftools.common.lims;

import static org.junit.Assert.assertEquals;

import org.junit.Test;

public class LimsGermlineReportingChoiceTest {

    @Test
    public void canExtractGermlineChoice() {
        assertEquals(LimsGermlineReportingChoice.REPORT_WITH_NOTIFICATION,
                LimsGermlineReportingChoice.fromLimsGermlineReportingChoiceString("1: Behandelbare toevalsbevindingen", "WIDE02991111T"));
        assertEquals(LimsGermlineReportingChoice.REPORT_WITH_NOTIFICATION,
                LimsGermlineReportingChoice.fromLimsGermlineReportingChoiceString("2: Alle toevalsbevindingen", "WIDE02991111T"));
        assertEquals(LimsGermlineReportingChoice.REPORT_WITHOUT_NOTIFICATION,
                LimsGermlineReportingChoice.fromLimsGermlineReportingChoiceString(
                        "3: Geen toevalsbevindingen; familie mag deze wel opvragen",
                        "WIDE02991111T"));
        assertEquals(LimsGermlineReportingChoice.REPORT_WITHOUT_NOTIFICATION,
                LimsGermlineReportingChoice.fromLimsGermlineReportingChoiceString("3: Geen toevalsbevindingen", "WIDE02991111T"));
        assertEquals(LimsGermlineReportingChoice.REPORT_WITHOUT_NOTIFICATION,
                LimsGermlineReportingChoice.fromLimsGermlineReportingChoiceString(
                        "4: Geen toevalsbevindingen; familie mag deze niet opvragen",
                        "WIDE02991111T"));
        assertEquals(LimsGermlineReportingChoice.NO_REPORTING,
                LimsGermlineReportingChoice.fromLimsGermlineReportingChoiceString("", "CPCT02991111T"));
        assertEquals(LimsGermlineReportingChoice.NO_REPORTING,
                LimsGermlineReportingChoice.fromLimsGermlineReportingChoiceString("", "DRUP02991111T"));
        assertEquals(LimsGermlineReportingChoice.NO_REPORTING,
                LimsGermlineReportingChoice.fromLimsGermlineReportingChoiceString("", "COLO02991111T"));
        assertEquals(LimsGermlineReportingChoice.REPORT_WITHOUT_NOTIFICATION,
                LimsGermlineReportingChoice.fromLimsGermlineReportingChoiceString("", "CORE02991111T"));
    }

    @Test(expected = IllegalStateException.class)
    public void hasUnknownGermlineChoice() {
        LimsGermlineReportingChoice.fromLimsGermlineReportingChoiceString("ALL", "WIDE02991111T");
    }
}