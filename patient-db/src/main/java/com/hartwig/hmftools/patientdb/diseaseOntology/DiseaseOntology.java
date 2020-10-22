package com.hartwig.hmftools.patientdb.diseaseOntology;

import java.io.File;
import java.util.logging.Level;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.jetbrains.annotations.NotNull;
import org.semanticweb.HermiT.Reasoner;
import org.semanticweb.owlapi.apibinding.OWLManager;
import org.semanticweb.owlapi.io.FileDocumentSource;
import org.semanticweb.owlapi.model.IRI;
import org.semanticweb.owlapi.model.MissingImportHandlingStrategy;
import org.semanticweb.owlapi.model.OWLClass;
import org.semanticweb.owlapi.model.OWLOntology;
import org.semanticweb.owlapi.model.OWLOntologyCreationException;
import org.semanticweb.owlapi.model.OWLOntologyLoaderConfiguration;
import org.semanticweb.owlapi.model.OWLOntologyManager;
import org.semanticweb.owlapi.reasoner.InferenceType;
import org.semanticweb.owlapi.reasoner.OWLReasoner;

public class DiseaseOntology {

    private static final Logger LOGGER = LogManager.getLogger(DiseaseOntology.class);

    private static final IRI DISEASE_IRI = IRI.create("http://purl.obolibrary.org/obo/DOID_4");

    public static void readDoid(@NotNull String doidFile) throws OWLOntologyCreationException {
        OWLOntology ontology = createOntology(doidFile);
        OWLReasoner reasoner = createReasoner(ontology);
        OWLClass a = createDiseaseMapping(ontology, reasoner);
    }

    private static OWLOntology createOntology(@NotNull String doidFile) throws OWLOntologyCreationException {
        java.util.logging.Logger logger = java.util.logging.Logger.getLogger("org.obolibrary.oboformat.parser.OBOFormatParser");
        logger.setLevel(Level.SEVERE);
        OWLOntologyManager ontologyManager = OWLManager.createOWLOntologyManager();
        OWLOntologyLoaderConfiguration config = new OWLOntologyLoaderConfiguration().setFollowRedirects(true)
                .setMissingImportHandlingStrategy(MissingImportHandlingStrategy.SILENT);

        return ontologyManager.loadOntologyFromOntologyDocument(new FileDocumentSource(new File(doidFile)), config);
    }

    private static OWLReasoner createReasoner(@NotNull OWLOntology ontology) {
        OWLReasoner reasoner = new Reasoner.ReasonerFactory().createReasoner(ontology);
        reasoner.precomputeInferences(InferenceType.CLASS_HIERARCHY);
        return reasoner;
    }

    private static OWLClass createDiseaseMapping(@NotNull OWLOntology ontology, @NotNull OWLReasoner reasoner) {

        OWLClass diseaseClass = ontology.getOWLOntologyManager().getOWLDataFactory().getOWLClass(DISEASE_IRI);
        LOGGER.info(diseaseClass);
        return diseaseClass;
        //        val diseases = setOf(diseaseClass) + subClasses(diseaseClass, reasoner)
        //        val synonyms = diseases.flatMap { disease -> getSynonyms(disease, ontology).map { Pair(it, disease) } }
        //        return diseases.associateBy { getLabel(it, ontology) } + synonyms.toMap()
    }
}
