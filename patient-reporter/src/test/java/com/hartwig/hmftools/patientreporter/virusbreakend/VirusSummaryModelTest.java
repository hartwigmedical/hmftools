package com.hartwig.hmftools.patientreporter.virusbreakend;

import static org.junit.Assert.*;

import java.util.Map;

import com.google.common.collect.Maps;
import com.google.common.collect.Sets;

import org.junit.Test;

public class VirusSummaryModelTest {

    @Test
    public void matchVirusToId() {
        Map<Integer, String> VirusIdMap = Maps.newHashMap();
        VirusIdMap.put(1, "virus1");
        VirusSummaryModel virusSummaryModel = new VirusSummaryModel(VirusIdMap);

        assertEquals("virus1", virusSummaryModel.findVirusSummary(1));
        assertNotEquals("virus2", virusSummaryModel.findVirusSummary(1));

        assertTrue(virusSummaryModel.mapIdtoVirusName(1));
        assertFalse(virusSummaryModel.mapIdtoVirusName(3));

    }
}