package com.hartwig.hmftools.knowledgebaseimporter.diseaseOntology

import org.semanticweb.HermiT.Reasoner
import org.semanticweb.owlapi.apibinding.OWLManager
import org.semanticweb.owlapi.io.FileDocumentSource
import org.semanticweb.owlapi.model.*
import org.semanticweb.owlapi.reasoner.InferenceType
import org.semanticweb.owlapi.reasoner.OWLReasoner
import java.io.File

class DiseaseOntology(fileLocation: String) {
    companion object {
        private val cancerIRI = IRI.create("http://purl.obolibrary.org/obo/DOID_162")

        private fun createOntology(fileLocation: String): OWLOntology {
            val ontologyManager = OWLManager.createOWLOntologyManager()
            val config = OWLOntologyLoaderConfiguration().setFollowRedirects(true).setMissingImportHandlingStrategy(
                    MissingImportHandlingStrategy.SILENT)
            return ontologyManager.loadOntologyFromOntologyDocument(FileDocumentSource(File(fileLocation)), config)
        }

        private fun createReasoner(ontology: OWLOntology): OWLReasoner {
            val reasoner = Reasoner.ReasonerFactory().createReasoner(ontology)
            reasoner.precomputeInferences(InferenceType.CLASS_HIERARCHY)
            return reasoner
        }

        private fun createCancerClassMapping(ontology: OWLOntology, reasoner: OWLReasoner): Map<String, OWLClass> {
            val cancerClass = ontology.owlOntologyManager.owlDataFactory.getOWLClass(cancerIRI)
            val cancerClasses = setOf(cancerClass) + reasoner.getSubClasses(cancerClass, false).flattened.filterNot { it.isOWLNothing }
            return cancerClasses.associateBy { getLabel(it, ontology) }
        }

        private fun getLabel(owlClass: OWLClass, ontology: OWLOntology): String {
            val dataFactory = ontology.owlOntologyManager.owlDataFactory
            return owlClass.getAnnotations(ontology, dataFactory.rdfsLabel).map { it.value as OWLLiteral }.map { it.literal }.first()
        }
    }

    private val ontology = createOntology(fileLocation)
    private val reasoner = createReasoner(ontology)
    val cancerClassMap = createCancerClassMapping(ontology, reasoner)
}
