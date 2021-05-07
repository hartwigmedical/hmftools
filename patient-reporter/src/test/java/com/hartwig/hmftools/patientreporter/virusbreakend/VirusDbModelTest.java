package com.hartwig.hmftools.patientreporter.virusbreakend;

import static org.junit.Assert.*;

import java.util.Map;

import com.google.common.collect.Maps;

import org.junit.Test;

public class VirusDbModelTest {

    @Test
    public void matchVirusToId() {
        Map<Integer, String> VirusIdMap = Maps.newHashMap();
        VirusIdMap.put(1, "virus1");
        VirusDbModel virusDbModel = new VirusDbModel(VirusIdMap);

        assertEquals("virus1", virusDbModel.findVirus(1));
        assertNotEquals("virus2", virusDbModel.findVirus(1));

        assertTrue(virusDbModel.mapIdtoVirusName(1));
        assertFalse(virusDbModel.mapIdtoVirusName(3));

    }
}