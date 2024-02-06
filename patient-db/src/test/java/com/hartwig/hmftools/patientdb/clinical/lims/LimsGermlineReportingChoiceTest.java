package com.hartwig.hmftools.patientdb.clinical.lims;

import static org.junit.Assert.assertEquals;

import com.hartwig.hmftools.patientdb.clinical.lims.cohort.LimsCohortConfig;
import com.hartwig.hmftools.patientdb.clinical.lims.cohort.LimsCohortTestFactory;

import org.junit.Test;

public class LimsGermlineReportingChoiceTest
{
    @Test
    public void canExtractGermlineLevel()
    {
        LimsCohortConfig cohortConfigCohort01 = LimsCohortTestFactory.createConfigForGermlineReporting("COHORT_01", true, true);
        assertEquals(LimsGermlineReportingLevel.REPORT_WITH_NOTIFICATION,
                LimsGermlineReportingLevel.fromLimsInputs(true, "1: Yes", "CHT01991111T", cohortConfigCohort01));
        assertEquals(LimsGermlineReportingLevel.REPORT_WITHOUT_NOTIFICATION,
                LimsGermlineReportingLevel.fromLimsInputs(true, "2: No", "CHT01991111T", cohortConfigCohort01));

        LimsCohortConfig cohortConfigCohort02 = LimsCohortTestFactory.createConfigForGermlineReporting("COHORT_02", true, true);
        assertEquals(LimsGermlineReportingLevel.REPORT_WITH_NOTIFICATION,
                LimsGermlineReportingLevel.fromLimsInputs(true, "1: Behandelbare toevalsbevindingen", "CHT022991111T", cohortConfigCohort02));
        assertEquals(LimsGermlineReportingLevel.REPORT_WITH_NOTIFICATION,
                LimsGermlineReportingLevel.fromLimsInputs(true, "2: Alle toevalsbevindingen", "CHT0202991111T", cohortConfigCohort02));
        assertEquals(LimsGermlineReportingLevel.REPORT_WITHOUT_NOTIFICATION,
                LimsGermlineReportingLevel.fromLimsInputs(true,
                        "3: Geen toevalsbevindingen; familie mag deze wel opvragen",
                        "CHT0202991111T",
                        cohortConfigCohort02));
        assertEquals(LimsGermlineReportingLevel.REPORT_WITHOUT_NOTIFICATION,
                LimsGermlineReportingLevel.fromLimsInputs(true, "3: Geen toevalsbevindingen", "CHT0202991111T", cohortConfigCohort02));
        assertEquals(LimsGermlineReportingLevel.REPORT_WITHOUT_NOTIFICATION,
                LimsGermlineReportingLevel.fromLimsInputs(true,
                        "4: Geen toevalsbevindingen; familie mag deze niet opvragen",
                        "CHT0202991111T",
                        cohortConfigCohort02));

        LimsCohortConfig cohortConfigCohort03 = LimsCohortTestFactory.createConfigForGermlineReporting("COHORT_03", false, false);
        assertEquals(LimsGermlineReportingLevel.NO_REPORTING,
                LimsGermlineReportingLevel.fromLimsInputs(false, "", "CHT0302991111T", cohortConfigCohort03));
        assertEquals(LimsGermlineReportingLevel.NO_REPORTING,
                LimsGermlineReportingLevel.fromLimsInputs(false, "", "CHT0302991111T", cohortConfigCohort03));

        LimsCohortConfig cohortConfigCohort04 = LimsCohortTestFactory.createConfigForGermlineReporting("COHORT_04", false, false);
        assertEquals(LimsGermlineReportingLevel.NO_REPORTING,
                LimsGermlineReportingLevel.fromLimsInputs(false, "", "CHT0402991111T", cohortConfigCohort04));

        LimsCohortConfig cohortConfigCohort05 = LimsCohortTestFactory.createConfigForGermlineReporting("COHORT_05", true, false);
        assertEquals(LimsGermlineReportingLevel.NO_REPORTING,
                LimsGermlineReportingLevel.fromLimsInputs(false, "", "CHT0502991111T", cohortConfigCohort05));
        assertEquals(LimsGermlineReportingLevel.NO_REPORTING,
                LimsGermlineReportingLevel.fromLimsInputs(false, "3: Geen toevalsbevindingen", "CHT0502991111T", cohortConfigCohort05));
    }

    @Test(expected = IllegalStateException.class)
    public void hasUnknownGermlineChoiceCohort05()
    {
        LimsCohortConfig cohortConfig = LimsCohortTestFactory.createConfigForGermlineReporting("CHT05", true, true);
        LimsGermlineReportingLevel.fromLimsInputs(true, "ALL", "CHT05991111T", cohortConfig);
    }

    @Test(expected = IllegalStateException.class)
    public void hasUnknownGermlineChoiceCohort06()
    {
        LimsCohortConfig cohortConfig = LimsCohortTestFactory.createConfigForGermlineReporting("CHT06", true, true);
        LimsGermlineReportingLevel.fromLimsInputs(true, "ALL", "CHT0602991111T", cohortConfig);
    }
}