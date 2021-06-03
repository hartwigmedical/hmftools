package com.hartwig.hmftools.patientreporter.virusbreakend;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertFalse;
import static org.junit.Assert.assertNotEquals;
import static org.junit.Assert.assertTrue;

import java.util.Map;

import com.google.common.collect.Maps;
import com.hartwig.hmftools.common.virus.VirusInterpretation;

import org.junit.Test;

public class VirusInterpretationModelTest {

    @Test
    public void canInterpretVirus() {
        Map<Integer, VirusInterpretation> speciesToInterpretationMap = Maps.newHashMap();
        speciesToInterpretationMap.put(1, VirusInterpretation.EBV);
        VirusInterpretationModel virusInterpretationModel = new VirusInterpretationModel(speciesToInterpretationMap);

        assertEquals(VirusInterpretation.EBV, virusInterpretationModel.interpretVirusSpecies(1));
        assertNotEquals(VirusInterpretation.HPV, virusInterpretationModel.interpretVirusSpecies(1));

        assertTrue(virusInterpretationModel.hasInterpretation(1));
        assertFalse(virusInterpretationModel.hasInterpretation(3));
    }
}