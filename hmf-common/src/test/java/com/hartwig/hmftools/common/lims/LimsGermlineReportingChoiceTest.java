package com.hartwig.hmftools.common.lims;

import static org.junit.Assert.assertEquals;

import java.util.Map;

import com.google.common.collect.Maps;
import com.hartwig.hmftools.common.lims.cohort.ImmutableLimsCohortConfigData;
import com.hartwig.hmftools.common.lims.cohort.ImmutableLimsCohortModel;
import com.hartwig.hmftools.common.lims.cohort.LimsCohortConfigData;
import com.hartwig.hmftools.common.lims.cohort.LimsCohortModel;

import org.jetbrains.annotations.NotNull;
import org.junit.Ignore;
import org.junit.Test;

public class LimsGermlineReportingChoiceTest {

    @NotNull
    private static LimsCohortModel buildTestCohortModel(@NotNull String cohortString) {
        Map<String, LimsCohortConfigData> cohortData = Maps.newHashMap();
        LimsCohortConfigData config = ImmutableLimsCohortConfigData.builder()
                .cohortId(cohortString)
                .hospitalId(true)
                .reportGermline(false)
                .reportGermlineFlag(false)
                .reportConclusion(false)
                .reportViral(false)
                .requireHospitalId(false)
                .requireHospitalPAId(false)
                .hospitalPersonsStudy(true)
                .hospitalPersonsRequester(false)
                .outputFile(false)
                .submission(false)
                .sidePanelInfo(false)
                .build();
        cohortData.put(cohortString, config);
        return ImmutableLimsCohortModel.builder().limsCohortMap(cohortData).build();
    }

    @Test
    @Ignore
    public void canExtractGermlineLevel() {
        LimsCohortModel cohortConfig = buildTestCohortModel("CPCT");
        assertEquals(LimsGermlineReportingLevel.REPORT_WITH_NOTIFICATION,
                LimsGermlineReportingLevel.fromLimsInputs(true, "1: Yes", "COREDB991111T", cohortConfig.queryCohortData("COREDB")));
        assertEquals(LimsGermlineReportingLevel.REPORT_WITHOUT_NOTIFICATION,
                LimsGermlineReportingLevel.fromLimsInputs(true, "2: No", "COREDB991111T", cohortConfig.queryCohortData("COREDB")));
        assertEquals(LimsGermlineReportingLevel.REPORT_WITH_NOTIFICATION,
                LimsGermlineReportingLevel.fromLimsInputs(true,
                        "1: Behandelbare toevalsbevindingen",
                        "WIDE02991111T",
                        cohortConfig.queryCohortData("WIDE")));
        assertEquals(LimsGermlineReportingLevel.REPORT_WITH_NOTIFICATION,
                LimsGermlineReportingLevel.fromLimsInputs(true,
                        "2: Alle toevalsbevindingen",
                        "WIDE02991111T",
                        cohortConfig.queryCohortData("WIDE")));
        assertEquals(LimsGermlineReportingLevel.REPORT_WITHOUT_NOTIFICATION,
                LimsGermlineReportingLevel.fromLimsInputs(true,
                        "3: Geen toevalsbevindingen; familie mag deze wel opvragen",
                        "WIDE02991111T",
                        cohortConfig.queryCohortData("WIDE")));
        assertEquals(LimsGermlineReportingLevel.REPORT_WITHOUT_NOTIFICATION,
                LimsGermlineReportingLevel.fromLimsInputs(true,
                        "3: Geen toevalsbevindingen",
                        "WIDE02991111T",
                        cohortConfig.queryCohortData("WIDE")));
        assertEquals(LimsGermlineReportingLevel.REPORT_WITHOUT_NOTIFICATION,
                LimsGermlineReportingLevel.fromLimsInputs(true,
                        "4: Geen toevalsbevindingen; familie mag deze niet opvragen",
                        "WIDE02991111T",
                        cohortConfig.queryCohortData("WIDE")));
        assertEquals(LimsGermlineReportingLevel.NO_REPORTING,
                LimsGermlineReportingLevel.fromLimsInputs(false, "", "CPCT02991111T", cohortConfig.queryCohortData("CPCT")));
        assertEquals(LimsGermlineReportingLevel.NO_REPORTING,
                LimsGermlineReportingLevel.fromLimsInputs(false, "", "DRUP02991111T", cohortConfig.queryCohortData("DRUP")));
        assertEquals(LimsGermlineReportingLevel.NO_REPORTING,
                LimsGermlineReportingLevel.fromLimsInputs(false, "", "COLO02991111T", cohortConfig.queryCohortData("CPCT")));
        assertEquals(LimsGermlineReportingLevel.REPORT_WITHOUT_NOTIFICATION,
                LimsGermlineReportingLevel.fromLimsInputs(true, "", "CORE02991111T", cohortConfig.queryCohortData("CORE")));

        assertEquals(LimsGermlineReportingLevel.NO_REPORTING,
                LimsGermlineReportingLevel.fromLimsInputs(false, "", "CORE02991111T", cohortConfig.queryCohortData("CORE")));
        assertEquals(LimsGermlineReportingLevel.NO_REPORTING,
                LimsGermlineReportingLevel.fromLimsInputs(false,
                        "3: Geen toevalsbevindingen",
                        "WIDE02991111T",
                        cohortConfig.queryCohortData("WIDE")));
    }

    @Ignore
    @Test(expected = IllegalStateException.class)
    public void hasUnknownGermlineChoice() {
        LimsCohortModel cohortConfig = buildTestCohortModel("CPCT");

        LimsGermlineReportingLevel.fromLimsInputs(true, "ALL", "WIDE02991111T", cohortConfig.queryCohortData("WIDE"));
        LimsGermlineReportingLevel.fromLimsInputs(true, "ALL", "COREDB991111T", cohortConfig.queryCohortData("COREDB"));
    }
}