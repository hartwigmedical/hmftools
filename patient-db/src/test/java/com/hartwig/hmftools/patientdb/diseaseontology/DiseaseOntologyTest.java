package com.hartwig.hmftools.patientdb.diseaseontology;

import java.io.FileNotFoundException;
import java.io.IOException;

import com.google.common.io.Resources;

import org.junit.Test;
import org.semanticweb.owlapi.model.OWLOntologyCreationException;

public class DiseaseOntologyTest {

    public static final String DOID_FILE_OWL = Resources.getResource("disease_ontology/doid.owl").getPath();
    public static final String DOID_FIL_JSON = Resources.getResource("disease_ontology/doid.json").getPath();

    @Test
    public void canLoadTest() throws OWLOntologyCreationException {
        try {
            DiseaseOntology.readDoid(DOID_FILE_OWL);
        } catch (OWLOntologyCreationException e) {
            throw new OWLOntologyCreationException("Could not load doid file!");
        }
    }


    @Test
    public void canLoadDoidJsonFile() throws IOException {
        try {
            DiseaseOntology.readDoidJsonFile(DOID_FIL_JSON);
        } catch (FileNotFoundException e) {
            throw new FileNotFoundException("Could not load doid file!");
        }
    }
}