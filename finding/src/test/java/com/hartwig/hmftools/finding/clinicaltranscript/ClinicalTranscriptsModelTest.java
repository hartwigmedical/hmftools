package com.hartwig.hmftools.finding.clinicaltranscript;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertNull;

import java.util.HashMap;
import java.util.Map;

import org.junit.Test;

public class ClinicalTranscriptsModelTest
{

    @Test
    public void canExtractClinicalTranscript()
    {
        Map<String, String> clinicalTranscriptMap = new HashMap<>();
        clinicalTranscriptMap.put("BRCA2", "NM_345");
        ClinicalTranscriptsModel clinicalTranscriptsModel = new ClinicalTranscriptsModel(clinicalTranscriptMap);

        assertEquals("NM_345", clinicalTranscriptsModel.findCanonicalTranscriptForGene("BRCA2"));
        assertNull(clinicalTranscriptsModel.findCanonicalTranscriptForGene("BRCA1"));
    }
}