package com.hartwig.hmftools.common.lims;

import static org.junit.Assert.assertEquals;

import org.junit.Test;

public class LimsGermlineReportingChoiceTest {

    @Test
    public void canExtractGermlineLevel() {
        assertEquals(LimsGermlineReportingLevel.REPORT_WITH_NOTIFICATION,
                LimsGermlineReportingLevel.fromLimsInputs(true, "1: Behandelbare toevalsbevindingen", "WIDE02991111T", LimsCohort.WIDE));
        assertEquals(LimsGermlineReportingLevel.REPORT_WITH_NOTIFICATION,
                LimsGermlineReportingLevel.fromLimsInputs(true, "2: Alle toevalsbevindingen", "WIDE02991111T", LimsCohort.WIDE));
        assertEquals(LimsGermlineReportingLevel.REPORT_WITHOUT_NOTIFICATION,
                LimsGermlineReportingLevel.fromLimsInputs(true,
                        "3: Geen toevalsbevindingen; familie mag deze wel opvragen",
                        "WIDE02991111T",
                        LimsCohort.WIDE));
        assertEquals(LimsGermlineReportingLevel.REPORT_WITHOUT_NOTIFICATION,
                LimsGermlineReportingLevel.fromLimsInputs(true, "3: Geen toevalsbevindingen", "WIDE02991111T", LimsCohort.WIDE));
        assertEquals(LimsGermlineReportingLevel.REPORT_WITHOUT_NOTIFICATION,
                LimsGermlineReportingLevel.fromLimsInputs(true,
                        "4: Geen toevalsbevindingen; familie mag deze niet opvragen",
                        "WIDE02991111T",
                        LimsCohort.WIDE));
        assertEquals(LimsGermlineReportingLevel.NO_REPORTING,
                LimsGermlineReportingLevel.fromLimsInputs(false, "", "CPCT02991111T", LimsCohort.CPCT));
        assertEquals(LimsGermlineReportingLevel.NO_REPORTING,
                LimsGermlineReportingLevel.fromLimsInputs(false, "", "DRUP02991111T", LimsCohort.DRUP));
        assertEquals(LimsGermlineReportingLevel.NO_REPORTING,
                LimsGermlineReportingLevel.fromLimsInputs(false, "", "COLO02991111T", LimsCohort.CPCT));
        assertEquals(LimsGermlineReportingLevel.REPORT_WITHOUT_NOTIFICATION,
                LimsGermlineReportingLevel.fromLimsInputs(true, "", "CORE02991111T", LimsCohort.CORE));

        assertEquals(LimsGermlineReportingLevel.NO_REPORTING,
                LimsGermlineReportingLevel.fromLimsInputs(false, "", "CORE02991111T", LimsCohort.CORE));
        assertEquals(LimsGermlineReportingLevel.NO_REPORTING,
                LimsGermlineReportingLevel.fromLimsInputs(false, "3: Geen toevalsbevindingen", "WIDE02991111T", LimsCohort.WIDE));
    }

    @Test(expected = IllegalStateException.class)
    public void hasUnknownGermlineChoice() {
        LimsGermlineReportingLevel.fromLimsInputs(true, "ALL", "WIDE02991111T", LimsCohort.WIDE);
    }
}