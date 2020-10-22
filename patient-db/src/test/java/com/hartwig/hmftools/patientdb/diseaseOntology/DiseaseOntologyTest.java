package com.hartwig.hmftools.patientdb.diseaseOntology;


import com.google.common.io.Resources;

import org.semanticweb.owlapi.model.OWLOntologyCreationException;

public class DiseaseOntologyTest {

    public static final String DOID_FILE = Resources.getResource("diseaseOntology/doid.owl").getPath();

    public static void canLoadDoidFile() throws OWLOntologyCreationException {
        try {
            DiseaseOntology.readDoid(DOID_FILE);
        } catch (OWLOntologyCreationException e) {
            throw new OWLOntologyCreationException("Could not load doid file!");
        }
    }

}