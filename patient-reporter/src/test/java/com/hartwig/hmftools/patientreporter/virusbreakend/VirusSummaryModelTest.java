package com.hartwig.hmftools.patientreporter.virusbreakend;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertFalse;
import static org.junit.Assert.assertNotEquals;
import static org.junit.Assert.assertTrue;

import java.util.Map;

import com.google.common.collect.Maps;

import org.junit.Test;

public class VirusSummaryModelTest {

    @Test
    public void matchVirusToId() {
        Map<Integer, String> VirusIdMap = Maps.newHashMap();
        VirusIdMap.put(1, "virus1");
        VirusSummaryModel virusSummaryModel = new VirusSummaryModel(VirusIdMap);

        assertEquals("virus1", virusSummaryModel.findVirusSummary(1));
        assertNotEquals("virus2", virusSummaryModel.findVirusSummary(1));

        assertTrue(virusSummaryModel.mapIdToVirusName(1));
        assertFalse(virusSummaryModel.mapIdToVirusName(3));
    }
}