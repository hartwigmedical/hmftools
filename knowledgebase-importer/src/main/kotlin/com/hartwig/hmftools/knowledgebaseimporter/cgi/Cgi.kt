package com.hartwig.hmftools.knowledgebaseimporter.cgi

import com.hartwig.hmftools.knowledgebaseimporter.Knowledgebase
import com.hartwig.hmftools.knowledgebaseimporter.diseaseOntology.DiseaseOntology
import com.hartwig.hmftools.knowledgebaseimporter.knowledgebases.KnowledgebaseSource
import com.hartwig.hmftools.knowledgebaseimporter.knowledgebases.RecordAnalyzer
import com.hartwig.hmftools.knowledgebaseimporter.output.*
import com.hartwig.hmftools.knowledgebaseimporter.readTSVRecords

class Cgi(variantsLocation: String, biomarkersLocation: String, diseaseOntology: DiseaseOntology,
          private val recordAnalyzer: RecordAnalyzer, treatmentTypeMap: Map<String, String>) :
        Knowledgebase, KnowledgebaseSource<CgiKnownKbRecord, CgiActionableRecord> {

    override val source = "cgi"
    override val knownVariants by lazy { recordAnalyzer.knownVariants(listOf(this)).distinct() }
    override val knownFusionPairs: List<FusionPair> by lazy { actionableKbRecords.flatMap { it.events }.filterIsInstance<FusionPair>().distinct() }
    override val promiscuousGenes: List<PromiscuousGene> by lazy { actionableKbRecords.flatMap { it.events }.filterIsInstance<PromiscuousGene>().distinct() }
    override val actionableVariants: List<ActionableVariantOutput> by lazy { actionableKbItems.filterIsInstance<ActionableVariantOutput>() }
    override val actionableCNVs: List<ActionableCNVOutput> by lazy { actionableKbItems.filterIsInstance<ActionableCNVOutput>() }
    override val actionableFusionPairs by lazy { actionableKbItems.filterIsInstance<ActionableFusionPairOutput>() }
    override val actionablePromiscuousGenes by lazy { actionableKbItems.filterIsInstance<ActionablePromiscuousGeneOutput>() }
    override val actionableRanges by lazy { actionableKbItems.filterIsInstance<ActionableGenomicRangeOutput>() }
    override val actionableGenes by lazy { actionableKbItems.filterIsInstance<ActionableGeneOutput>() }
    override val cancerTypes by lazy {
        actionableKbRecords.flatMap { it.actionability }.map { it.cancerType }
                .associateBy({ it }, { diseaseOntology.findDoidsForCancerType(it) })
    }
    override val knownKbRecords by lazy { readTSVRecords(variantsLocation) { CgiKnownKbRecord(it) }.filterNotNull() }
    override val actionableKbRecords by lazy {
        readTSVRecords(biomarkersLocation) { CgiActionableRecord(it, treatmentTypeMap) }
    }
    private val actionableKbItems by lazy { recordAnalyzer.actionableItems(listOf(this)).distinct() }
}
