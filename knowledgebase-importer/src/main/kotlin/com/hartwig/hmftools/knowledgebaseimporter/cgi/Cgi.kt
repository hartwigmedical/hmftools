package com.hartwig.hmftools.knowledgebaseimporter.cgi

import com.hartwig.hmftools.knowledgebaseimporter.Knowledgebase
import com.hartwig.hmftools.knowledgebaseimporter.diseaseOntology.DiseaseOntology
import com.hartwig.hmftools.knowledgebaseimporter.knowledgebases.KnowledgebaseSource
import com.hartwig.hmftools.knowledgebaseimporter.knowledgebases.RecordAnalyzer
import com.hartwig.hmftools.knowledgebaseimporter.output.*
import com.hartwig.hmftools.knowledgebaseimporter.readTSVRecords
import htsjdk.samtools.reference.IndexedFastaSequenceFile

class Cgi(variantsLocation: String, biomarkersLocation: String, transvarLocation: String, diseaseOntology: DiseaseOntology,
          private val reference: IndexedFastaSequenceFile) :
        Knowledgebase, KnowledgebaseSource<CgiKnownKbRecord, CgiActionableRecord> {

    override val source = "cgi"
    override val knownVariants by lazy { RecordAnalyzer(transvarLocation, reference).knownVariants(listOf(this)).distinct() }
    override val knownFusionPairs: List<FusionPair> by lazy { actionableKbRecords.flatMap { it.events }.filterIsInstance<FusionPair>().distinct() }
    override val promiscuousGenes: List<PromiscuousGene> by lazy { actionableKbRecords.flatMap { it.events }.filterIsInstance<PromiscuousGene>().distinct() }
    override val actionableVariants: List<ActionableVariantOutput> by lazy { actionableKbVariants.map { it.toActionableOutput() }.filterIsInstance<ActionableVariantOutput>() }
    override val actionableCNVs: List<ActionableCNVOutput> by lazy { actionableKbVariants.map { it.toActionableOutput() }.filterIsInstance<ActionableCNVOutput>() }
    override val actionableFusions: List<ActionableFusionOutput> by lazy { actionableKbVariants.map { it.toActionableOutput() }.filterIsInstance<ActionableFusionOutput>() }
    override val cancerTypes by lazy {
        actionableKbRecords.flatMap { it.actionability }.map { it.cancerType }
                .associateBy({ it }, { diseaseOntology.findDoidsForCancerType(it) })
    }
    override val knownKbRecords by lazy { readTSVRecords(variantsLocation) { CgiKnownKbRecord(it) }.filterNotNull() }
    override val actionableKbRecords by lazy { readTSVRecords(biomarkersLocation) { CgiActionableRecord(it) }.filterNotNull() }
    val actionableKbVariants by lazy { RecordAnalyzer(transvarLocation, reference).actionableItems(listOf(this)).distinct() }
}
