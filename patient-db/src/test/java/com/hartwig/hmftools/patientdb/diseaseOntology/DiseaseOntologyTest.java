package com.hartwig.hmftools.patientdb.diseaseOntology;

import com.google.common.io.Resources;

import org.junit.Ignore;
import org.junit.Test;
import org.semanticweb.owlapi.model.OWLOntologyCreationException;

public class DiseaseOntologyTest {

    public static final String DOID_FILE = Resources.getResource("diseaseOntology/doid.owl").getPath();

    @Test
    public void canLoadTest() throws OWLOntologyCreationException {
        try {
            DiseaseOntology.readDoid(DOID_FILE);
        } catch (OWLOntologyCreationException e) {
            throw new OWLOntologyCreationException("Could not load doid file!");
        }
    }

}