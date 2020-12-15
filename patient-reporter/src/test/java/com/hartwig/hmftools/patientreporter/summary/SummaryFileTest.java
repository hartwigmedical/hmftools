package com.hartwig.hmftools.patientreporter.summary;

import static org.junit.Assert.assertEquals;

import java.io.IOException;
import java.util.Map;

import com.google.common.collect.Maps;
import com.google.common.io.Resources;
import com.hartwig.hmftools.common.lims.cohort.ImmutableLimsCohortConfigData;
import com.hartwig.hmftools.common.lims.cohort.ImmutableLimsCohortModel;
import com.hartwig.hmftools.common.lims.cohort.LimsCohortConfigData;
import com.hartwig.hmftools.common.lims.cohort.LimsCohortModel;

import org.jetbrains.annotations.NotNull;
import org.junit.Test;

public class SummaryFileTest {

    private static final String SAMPLE_SUMMARY_TSV = Resources.getResource("sample_summary/sample_summary.tsv").getPath();

    @Test
    public void summaryFromCSVWithNewLines() throws IOException {
        SummaryModel summaryModel = SummaryFile.buildFromTsv(SAMPLE_SUMMARY_TSV);
        assertEquals(1, summaryModel.summaryCount());

        String summary = summaryModel.findSummaryForSample("sample", buildTestCohortModel("CORE", true).queryCohortData("CORE", "CORE01010000T"));

        assertEquals(3, summary.split("\n").length);
    }

    @NotNull
    private static LimsCohortModel buildTestCohortModel(@NotNull String cohortString, boolean requireConclusion) {
        Map<String, LimsCohortConfigData> cohortData = Maps.newHashMap();
        LimsCohortConfigData config = ImmutableLimsCohortConfigData.builder()
                .cohortId(cohortString)
                .hospitalId(true)
                .reportGermline(false)
                .reportGermlineFlag(false)
                .reportConclusion(requireConclusion)
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
}