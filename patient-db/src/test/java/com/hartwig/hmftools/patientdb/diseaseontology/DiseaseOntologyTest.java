package com.hartwig.hmftools.patientdb.diseaseontology;

import static org.junit.Assert.assertEquals;

import java.io.FileNotFoundException;
import java.io.IOException;
import java.util.List;

import com.google.common.io.Resources;

import org.junit.Test;

public class DiseaseOntologyTest {

    private static final String DOID_FIL_JSON = Resources.getResource("disease_ontology/doid.json").getPath();

    @Test
    public void canLoadDoidJsonFile() throws IOException {
        try {
            List<Doid> doids = DiseaseOntology.readDoidJsonFile(DOID_FIL_JSON);
            assertEquals(2, doids.size());

        } catch (FileNotFoundException e) {
            throw new FileNotFoundException("Could not load doid file!");
        }
    }
}