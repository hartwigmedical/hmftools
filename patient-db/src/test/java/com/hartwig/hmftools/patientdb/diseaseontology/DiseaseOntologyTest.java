package com.hartwig.hmftools.patientdb.diseaseontology;

import java.io.FileNotFoundException;
import java.io.IOException;

import com.google.common.io.Resources;

import org.junit.Test;

public class DiseaseOntologyTest {

    private static final String DOID_FIL_JSON = Resources.getResource("disease_ontology/doid.json").getPath();

    @Test
    public void canLoadDoidJsonFile() throws IOException {
        try {
            DiseaseOntology.readDoidJsonFile(DOID_FIL_JSON);
        } catch (FileNotFoundException e) {
            throw new FileNotFoundException("Could not load doid file!");
        }
    }
}