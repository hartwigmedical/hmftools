package com.hartwig.hmftools.patientdb.clinical.lims;

import static org.junit.Assert.assertEquals;

import com.hartwig.hmftools.patientdb.clinical.lims.cohort.LimsCohortConfig;
import com.hartwig.hmftools.patientdb.clinical.lims.cohort.LimsCohortTestFactory;

import org.junit.Test;

public class LimsGermlineReportingChoiceTest {

    @Test
    public void canExtractGermlineLevel() {
        LimsCohortConfig cohortConfigCOREDB = LimsCohortTestFactory.createConfigForGermlineReporting("COREDB", true, true);
        assertEquals(LimsGermlineReportingLevel.REPORT_WITH_NOTIFICATION,
                LimsGermlineReportingLevel.fromLimsInputs(true, "1: Yes", "COREDB991111T", cohortConfigCOREDB));
        assertEquals(LimsGermlineReportingLevel.REPORT_WITHOUT_NOTIFICATION,
                LimsGermlineReportingLevel.fromLimsInputs(true, "2: No", "COREDB991111T", cohortConfigCOREDB));

        LimsCohortConfig cohortConfigWIDE = LimsCohortTestFactory.createConfigForGermlineReporting("WIDE", true, true);
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

        LimsCohortConfig cohortConfigCPCT = LimsCohortTestFactory.createConfigForGermlineReporting("CPCT", false, false);
        assertEquals(LimsGermlineReportingLevel.NO_REPORTING,
                LimsGermlineReportingLevel.fromLimsInputs(false, "", "CPCT02991111T", cohortConfigCPCT));
        assertEquals(LimsGermlineReportingLevel.NO_REPORTING,
                LimsGermlineReportingLevel.fromLimsInputs(false, "", "COLO02991111T", cohortConfigCPCT));

        LimsCohortConfig cohortConfigDRUP = LimsCohortTestFactory.createConfigForGermlineReporting("DRUP", false, false);
        assertEquals(LimsGermlineReportingLevel.NO_REPORTING,
                LimsGermlineReportingLevel.fromLimsInputs(false, "", "DRUP02991111T", cohortConfigDRUP));

        LimsCohortConfig cohortConfigCORE = LimsCohortTestFactory.createConfigForGermlineReporting("CORE", true, false);
        assertEquals(LimsGermlineReportingLevel.NO_REPORTING,
                LimsGermlineReportingLevel.fromLimsInputs(false, "", "CORE02991111T", cohortConfigCORE));
        assertEquals(LimsGermlineReportingLevel.NO_REPORTING,
                LimsGermlineReportingLevel.fromLimsInputs(false, "3: Geen toevalsbevindingen", "WIDE02991111T", cohortConfigCORE));
    }

    @Test(expected = IllegalStateException.class)
    public void hasUnknownGermlineChoiceCOREDB() {
        LimsCohortConfig cohortConfigCOREDB = LimsCohortTestFactory.createConfigForGermlineReporting("COREDB", true, true);
        LimsGermlineReportingLevel.fromLimsInputs(true, "ALL", "COREDB991111T", cohortConfigCOREDB);
    }

    @Test(expected = IllegalStateException.class)
    public void hasUnknownGermlineChoiceWIDE() {
        LimsCohortConfig cohortConfigWIDE = LimsCohortTestFactory.createConfigForGermlineReporting("WIDE", true, true);
        LimsGermlineReportingLevel.fromLimsInputs(true, "ALL", "WIDE02991111T", cohortConfigWIDE);
    }
}