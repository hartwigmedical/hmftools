package com.hartwig.hmftools.common.lims;

import static org.junit.Assert.assertEquals;

import com.hartwig.hmftools.common.lims.cohort.LimsCohortConfig;

import org.junit.Test;

public class LimsGermlineReportingChoiceTest {

    @Test
    public void canExtractGermlineLevel() {
        LimsCohortConfig cohortConfigCOREDB =
                LimsTestUtil.createTestCohortConfig("COREDB", true, true, true, false, true, true, true, true, false, true, false, true);
        assertEquals(LimsGermlineReportingLevel.REPORT_WITH_NOTIFICATION,
                LimsGermlineReportingLevel.fromLimsInputs(true, "1: Yes", "COREDB991111T", cohortConfigCOREDB));
        assertEquals(LimsGermlineReportingLevel.REPORT_WITHOUT_NOTIFICATION,
                LimsGermlineReportingLevel.fromLimsInputs(true, "2: No", "COREDB991111T", cohortConfigCOREDB));

        LimsCohortConfig cohortConfigWIDE =
                LimsTestUtil.createTestCohortConfig("WIDE", true, true, true, true, true, false, true, true, false, false, false, true);
        assertEquals(LimsGermlineReportingLevel.REPORT_WITH_NOTIFICATION,
                LimsGermlineReportingLevel.fromLimsInputs(true, "1: Behandelbare toevalsbevindingen", "WIDE02991111T", cohortConfigWIDE));
        assertEquals(LimsGermlineReportingLevel.REPORT_WITH_NOTIFICATION,
                LimsGermlineReportingLevel.fromLimsInputs(true, "2: Alle toevalsbevindingen", "WIDE02991111T", cohortConfigWIDE));
        assertEquals(LimsGermlineReportingLevel.REPORT_WITHOUT_NOTIFICATION,
                LimsGermlineReportingLevel.fromLimsInputs(true,
                        "3: Geen toevalsbevindingen; familie mag deze wel opvragen",
                        "WIDE02991111T",
                        cohortConfigWIDE));
        assertEquals(LimsGermlineReportingLevel.REPORT_WITHOUT_NOTIFICATION,
                LimsGermlineReportingLevel.fromLimsInputs(true, "3: Geen toevalsbevindingen", "WIDE02991111T", cohortConfigWIDE));
        assertEquals(LimsGermlineReportingLevel.REPORT_WITHOUT_NOTIFICATION,
                LimsGermlineReportingLevel.fromLimsInputs(true,
                        "4: Geen toevalsbevindingen; familie mag deze niet opvragen",
                        "WIDE02991111T",
                        cohortConfigWIDE));

        LimsCohortConfig cohortConfigCPCT =
                LimsTestUtil.createTestCohortConfig("CPCT", true, false, false, false, false, false, false, true, false, false, false, false);
        assertEquals(LimsGermlineReportingLevel.NO_REPORTING,
                LimsGermlineReportingLevel.fromLimsInputs(false, "", "CPCT02991111T", cohortConfigCPCT));
        assertEquals(LimsGermlineReportingLevel.NO_REPORTING,
                LimsGermlineReportingLevel.fromLimsInputs(false, "", "COLO02991111T", cohortConfigCPCT));

        LimsCohortConfig cohortConfigDRUP =
                LimsTestUtil.createTestCohortConfig("DRUP", true, false, false, false, false, false, false, true, false, false, false, false);
        assertEquals(LimsGermlineReportingLevel.NO_REPORTING,
                LimsGermlineReportingLevel.fromLimsInputs(false, "", "DRUP02991111T", cohortConfigDRUP));

        LimsCohortConfig cohortConfigCORE =
                LimsTestUtil.createTestCohortConfig("CORE", true, true, false, true, true, true, true, false, true, true, true, true);
        assertEquals(LimsGermlineReportingLevel.NO_REPORTING,
                LimsGermlineReportingLevel.fromLimsInputs(false, "", "CORE02991111T", cohortConfigCORE));
        assertEquals(LimsGermlineReportingLevel.NO_REPORTING,
                LimsGermlineReportingLevel.fromLimsInputs(false, "3: Geen toevalsbevindingen", "WIDE02991111T", cohortConfigCORE));
    }

    @Test(expected = IllegalStateException.class)
    public void hasUnknownGermlineChoice() {
        LimsCohortConfig cohortConfigCOREDB =
                LimsTestUtil.createTestCohortConfig("COREDB", true, true, true, false, true, true, true, true, false, true, false, true);
        LimsCohortConfig cohortConfigWIDE =
                LimsTestUtil.createTestCohortConfig("WIDE", true, true, true, true, true, false, true, true, false, false, false, true);
        LimsCohortConfig cohortConfigCORE =
                LimsTestUtil.createTestCohortConfig("CORE", true, true, false, true, true, true, true, false, true, true, true, true);

        LimsGermlineReportingLevel.fromLimsInputs(true, "ALL", "WIDE02991111T", cohortConfigWIDE);
        LimsGermlineReportingLevel.fromLimsInputs(true, "ALL", "COREDB991111T", cohortConfigCOREDB);
        LimsGermlineReportingLevel.fromLimsInputs(true, "a", "CORE02991111T", cohortConfigCORE);
    }
}