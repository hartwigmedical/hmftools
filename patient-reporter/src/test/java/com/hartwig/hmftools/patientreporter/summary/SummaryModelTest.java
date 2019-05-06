package com.hartwig.hmftools.patientreporter.summary;

import static org.junit.Assert.*;

import java.util.Map;

import com.google.common.collect.Maps;

import org.junit.Test;

public class SummaryModelTest {

    @Test
    public void sampleArePresentInSummaryFile() {

        Map<String, String> summaryToSampleMap = Maps.newHashMap();
        summaryToSampleMap.put("CPCT11111111T", "this is a test summary");
        SummaryModel summaryModel = new SummaryModel(summaryToSampleMap);

        assertTrue(summaryModel.sampleIdPresentInSummaryFile("CPCT11111111T"));
        assertFalse(summaryModel.sampleIdPresentInSummaryFile("CPCT"));
    }

    @Test
    public void canExtractSummaryOfSample() {
        Map<String, String> summaryToSampleMap = Maps.newHashMap();
        summaryToSampleMap.put("CPCT11111111T", "this is a test summary");
        SummaryModel summaryModel = new SummaryModel(summaryToSampleMap);

        assertEquals("this is a test summary", summaryModel.extractSummarySampleId("CPCT11111111T"));
        assertNotEquals("this is a test summary", summaryModel.extractSummarySampleId("CPCT11111111"));
    }
}