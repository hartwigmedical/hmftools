package com.hartwig.hmftools.knowledgebaseimporter.civic

import com.google.common.collect.ArrayListMultimap
import com.google.common.collect.Multimap
import com.hartwig.hmftools.apiclients.civic.api.CivicApiWrapper
import com.hartwig.hmftools.knowledgebaseimporter.Knowledgebase
import com.hartwig.hmftools.knowledgebaseimporter.diseaseOntology.DiseaseOntology
import com.hartwig.hmftools.knowledgebaseimporter.knowledgebases.ActionableRecord
import com.hartwig.hmftools.knowledgebaseimporter.knowledgebases.KnowledgebaseSource
import com.hartwig.hmftools.knowledgebaseimporter.knowledgebases.RecordAnalyzer
import com.hartwig.hmftools.knowledgebaseimporter.output.*
import com.hartwig.hmftools.knowledgebaseimporter.readTSVRecords

class Civic(variantsLocation: String, evidenceLocation: String, diseaseOntology: DiseaseOntology,
            private val recordAnalyzer: RecordAnalyzer, treatmentTypeMap: Map<String, String>) :
        Knowledgebase, KnowledgebaseSource<CivicRecord, ActionableRecord> {

    override val source = "civic"
    override val knownVariants by lazy { recordAnalyzer.knownVariants(listOf(this)).distinct() }
    override val knownFusionPairs: List<FusionPair> by lazy { knownKbRecords.flatMap { it.events }.filterIsInstance<FusionPair>().distinct() }
    override val promiscuousGenes: List<PromiscuousGene> by lazy { knownKbRecords.flatMap { it.events }.filterIsInstance<PromiscuousGene>().distinct() }
    override val actionableVariants: List<ActionableVariantOutput> by lazy { actionableKbItems.filterIsInstance<ActionableVariantOutput>() }
    override val actionableCNVs: List<ActionableCNVOutput> by lazy { actionableKbItems.filterIsInstance<ActionableCNVOutput>() }
    override val actionableFusionPairs by lazy { actionableKbItems.filterIsInstance<ActionableFusionPairOutput>() }
    override val actionablePromiscuousGenes by lazy { actionableKbItems.filterIsInstance<ActionablePromiscuousGeneOutput>() }
    override val actionableRanges by lazy { actionableKbItems.filterIsInstance<ActionableGenomicRangeOutput>() }
    override val cancerTypes by lazy {
        actionableKbRecords.flatMap { it.cancerDoids.entries }
                .associateBy({ it.key }, { diseaseOntology.findDoidsForCancerType(it.key) + diseaseOntology.findDoidsForDoid(it.value) })
    }

    override val knownKbRecords by lazy { preProcessCivicRecords(variantsLocation, evidenceLocation, treatmentTypeMap) }
    override val actionableKbRecords by lazy { knownKbRecords }
    private val actionableKbItems by lazy { recordAnalyzer.actionableItems(listOf(this)).distinct() }

    private fun preProcessCivicRecords(variantFileLocation: String, evidenceFileLocation: String,
                                       treatmentTypeMap: Map<String, String>): List<CivicRecord> {
        val variantEvidenceMap = readEvidenceMap(evidenceFileLocation, treatmentTypeMap)
        return readTSVRecords(variantFileLocation) { CivicRecord(it, variantEvidenceMap) }
    }

    private fun readEvidenceMap(evidenceLocation: String, treatmentTypeMap: Map<String, String>): Multimap<String, CivicEvidence> {
        val civicApi = CivicApiWrapper()
        val drugInteractionMap = civicApi.drugInteractionMap
        val evidenceMap = ArrayListMultimap.create<String, CivicEvidence>()
        readTSVRecords(evidenceLocation) { csvRecord ->
            evidenceMap.put(csvRecord["variant_id"], CivicEvidence(csvRecord, drugInteractionMap, treatmentTypeMap))
        }
        civicApi.releaseResources()
        return evidenceMap
    }
}
