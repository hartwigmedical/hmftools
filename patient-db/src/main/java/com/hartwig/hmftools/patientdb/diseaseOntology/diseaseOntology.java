package com.hartwig.hmftools.patientdb.diseaseOntology;


import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.semanticweb.owlapi.model.IRI;

public class diseaseOntology {
    private static final Logger LOGGER = LogManager.getLogger(diseaseOntology.class);

    private static final IRI DISEASE_IRI = IRI.create("http://purl.obolibrary.org/obo/DOID_4");

}
