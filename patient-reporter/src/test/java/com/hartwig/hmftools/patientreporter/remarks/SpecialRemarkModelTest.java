package com.hartwig.hmftools.patientreporter.remarks;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertFalse;
import static org.junit.Assert.assertNotEquals;
import static org.junit.Assert.assertTrue;

import java.util.Map;

import com.google.common.collect.Maps;

import org.junit.Test;

public class SpecialRemarkModelTest {

    @Test
    public void sampleArePresentInSpecialRemarkModel() {
        Map<String, String> specialRemarkToSampleMap = Maps.newHashMap();
        specialRemarkToSampleMap.put("sample", "this is a test summary");
        SpecialRemarkModel specialRemarkModel = new SpecialRemarkModel(specialRemarkToSampleMap);

        assertTrue(specialRemarkModel.samplePresentInSpecialRemarks("sample"));
        assertFalse(specialRemarkModel.samplePresentInSpecialRemarks("sample2"));
    }

    @Test
    public void canExtractSpecialRemarkOfSample() {
        Map<String, String> specialRemarkToSampleMap = Maps.newHashMap();
        specialRemarkToSampleMap.put("sample", "this is a test summary");
        SpecialRemarkModel specialRemarkModel = new SpecialRemarkModel(specialRemarkToSampleMap);

        assertEquals("this is a test summary", specialRemarkModel.findSpecialRemarkForSample("sample"));
        assertNotEquals("this is a test summary", specialRemarkModel.findSpecialRemarkForSample("sample2"));
    }
}