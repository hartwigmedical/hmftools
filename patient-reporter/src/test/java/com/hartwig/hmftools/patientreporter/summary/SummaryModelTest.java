package com.hartwig.hmftools.patientreporter.summary;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertFalse;
import static org.junit.Assert.assertNotEquals;
import static org.junit.Assert.assertTrue;

import java.util.Map;

import com.google.common.collect.Maps;
import com.hartwig.hmftools.common.lims.cohort.LimsCohortConfig;
import com.hartwig.hmftools.common.lims.cohort.LimsCohortTestFactory;

import org.junit.Test;

public class SummaryModelTest {

    @Test
    public void sampleArePresentInSummaryModel() {
        Map<String, String> summaryToSampleMap = Maps.newHashMap();
        summaryToSampleMap.put("sample", "this is a test summary");
        SummaryModel summaryModel = new SummaryModel(summaryToSampleMap);

        assertTrue(summaryModel.samplePresentInSummaries("sample"));
        assertFalse(summaryModel.samplePresentInSummaries("sample2"));
    }

    @Test
    public void canExtractSummaryOfSample() {
        Map<String, String> summaryToSampleMap = Maps.newHashMap();
        summaryToSampleMap.put("sample", "this is a test summary");
        SummaryModel summaryModel = new SummaryModel(summaryToSampleMap);

        LimsCohortConfig cohortConfig = LimsCohortTestFactory.createCOLOCohortConfig();

        assertEquals("this is a test summary", summaryModel.findSummaryForSample("sample", cohortConfig));
        assertNotEquals("this is a test summary", summaryModel.findSummaryForSample("sample2", cohortConfig));
    }
}