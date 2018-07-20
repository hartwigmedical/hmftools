package com.hartwig.hmftools.knowledgebaseimporter.diseaseOntology

import org.semanticweb.HermiT.Reasoner
import org.semanticweb.owlapi.apibinding.OWLManager
import org.semanticweb.owlapi.io.FileDocumentSource
import org.semanticweb.owlapi.model.*
import org.semanticweb.owlapi.reasoner.InferenceType
import org.semanticweb.owlapi.reasoner.OWLReasoner
import java.io.File
import java.util.logging.Level
import java.util.logging.Logger

class DiseaseOntology(fileLocation: String) {
    companion object {
        private val diseaseIRI = IRI.create("http://purl.obolibrary.org/obo/DOID_4")

        private fun createOntology(fileLocation: String): OWLOntology {
            val logger = Logger.getLogger("org.obolibrary.oboformat.parser.OBOFormatParser")
            logger.level = Level.SEVERE
            val ontologyManager = OWLManager.createOWLOntologyManager()
            val config = OWLOntologyLoaderConfiguration().setFollowRedirects(true)
                    .setMissingImportHandlingStrategy(MissingImportHandlingStrategy.SILENT)
            return ontologyManager.loadOntologyFromOntologyDocument(FileDocumentSource(File(fileLocation)), config)
        }

        private fun createReasoner(ontology: OWLOntology): OWLReasoner {
            val reasoner = Reasoner.ReasonerFactory().createReasoner(ontology)
            reasoner.precomputeInferences(InferenceType.CLASS_HIERARCHY)
            return reasoner
        }

        private fun createDiseaseMapping(ontology: OWLOntology, reasoner: OWLReasoner): Map<String, OWLClass> {
            val diseaseClass = ontology.owlOntologyManager.owlDataFactory.getOWLClass(diseaseIRI)
            val diseases = setOf(diseaseClass) + subClasses(diseaseClass, reasoner)
            val synonyms = diseases.flatMap { disease -> getSynonyms(disease, ontology).map { Pair(it, disease) } }
            return diseases.associateBy { getLabel(it, ontology) } + synonyms.toMap()
        }

        private fun subClasses(owlClass: OWLClass, reasoner: OWLReasoner): Set<OWLClass> {
            return reasoner.getSubClasses(owlClass, false).flattened.filterNot { it.isOWLNothing || it.isOWLThing }.toSet()
        }

        private fun superClasses(owlClass: OWLClass, reasoner: OWLReasoner): Set<OWLClass> {
            return reasoner.getSuperClasses(owlClass, false).flattened.filterNot { it.isOWLNothing || it.isOWLThing }.toSet()
        }

        private fun getLabel(owlClass: OWLClass, ontology: OWLOntology): String {
            val dataFactory = ontology.owlOntologyManager.owlDataFactory
            return owlClass.getAnnotations(ontology, dataFactory.rdfsLabel)
                    .map { it.value as OWLLiteral }
                    .map { it.literal.toLowerCase().trim() }
                    .first()
        }

        private fun getSynonyms(owlClass: OWLClass, ontology: OWLOntology): Set<String> {
            return owlClass.getAnnotations(ontology)
                    .filter { it.property.toString().contains("hasExactSynonym") }
                    .map { it.value as OWLLiteral }
                    .map { it.literal.toLowerCase().trim() }
                    .toSet()
        }
    }

    private val ontology by lazy { createOntology(fileLocation) }
    private val reasoner by lazy { createReasoner(ontology) }
    private val diseaseNameToClass by lazy { createDiseaseMapping(ontology, reasoner) }

    fun isPresent(cancerType: String): Boolean {
        return diseaseNameToClass.containsKey(cancerType.toLowerCase().trim())
    }

    fun findDoidsForCancerType(cancerType: String): Set<String> {
        val cancerClass = diseaseNameToClass[cancerType.toLowerCase().trim()]
        cancerClass ?: return emptySet()
        return findDoids(cancerClass)
    }

    fun findDoidsForDoid(cancerDoid: String): Set<String> {
        if (cancerDoid.isEmpty()) return emptySet()
        val cancerIRI = IRI.create("http://purl.obolibrary.org/obo/DOID_$cancerDoid")
        return findDoids(ontology.owlOntologyManager.owlDataFactory.getOWLClass(cancerIRI))
    }

    private fun findDoids(cancerClass: OWLClass): Set<String> {
        val relevantClasses = setOf(cancerClass) + subClasses(cancerClass, reasoner) + superClasses(cancerClass, reasoner)
        return relevantClasses.map { it.toString().substringAfter("DOID_").substringBefore('>') }
                .filterNot { it.isBlank() }
                .toSet()
    }
}
