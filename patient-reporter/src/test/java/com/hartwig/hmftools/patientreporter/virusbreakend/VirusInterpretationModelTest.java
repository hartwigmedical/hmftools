package com.hartwig.hmftools.patientreporter.virusbreakend;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertFalse;
import static org.junit.Assert.assertNotEquals;
import static org.junit.Assert.assertTrue;

import java.util.Map;

import com.google.common.collect.Maps;

import org.junit.Test;

public class VirusInterpretationModelTest {

    @Test
    public void canInterpretVirus() {
        Map<Integer, String> speciesToInterpretationMap = Maps.newHashMap();
        speciesToInterpretationMap.put(1, "virus1");
        VirusInterpretationModel virusInterpretationModel = new VirusInterpretationModel(speciesToInterpretationMap);

        assertEquals("virus1", virusInterpretationModel.interpretVirusSpecies(1));
        assertNotEquals("virus2", virusInterpretationModel.interpretVirusSpecies(1));

        assertTrue(virusInterpretationModel.hasInterpretation(1));
        assertFalse(virusInterpretationModel.hasInterpretation(3));
    }
}