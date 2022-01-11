package com.hartwig.hmftools.virusinterpreter.algo;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertFalse;
import static org.junit.Assert.assertTrue;

import java.util.Map;

import com.google.common.collect.Maps;
import com.hartwig.hmftools.common.virus.VirusLikelihoodType;

import org.junit.Test;

public class VirusReportingDbModelTest {

    @Test
    public void canInterpretVirus() {
        Map<Integer, VirusReportingDb> speciesToInterpretationMap = Maps.newHashMap();

        VirusReportingDb virusWhitelist = ImmutableVirusReportingDb.builder()
                .virusInterpretation("EBV")
                .integratedMinimalCoverage(null)
                .nonIntegratedMinimalCoverage(null)
                .virusDriverLikelihoodType(VirusLikelihoodType.HIGH)
                .build();

        speciesToInterpretationMap.put(1, virusWhitelist);
        VirusReportingDbModel virusInterpretationModel = new VirusReportingDbModel(speciesToInterpretationMap);

        assertEquals("EBV", virusInterpretationModel.interpretVirusSpecies(1));

        assertTrue(virusInterpretationModel.hasInterpretation(1));
        assertFalse(virusInterpretationModel.hasInterpretation(3));
    }
}