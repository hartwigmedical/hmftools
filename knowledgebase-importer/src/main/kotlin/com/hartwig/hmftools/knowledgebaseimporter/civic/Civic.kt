package com.hartwig.hmftools.knowledgebaseimporter.civic

import com.google.common.collect.ArrayListMultimap
import com.google.common.collect.Multimap
import com.hartwig.hmftools.apiclients.civic.api.CivicApiWrapper
import com.hartwig.hmftools.extensions.csv.CsvReader
import com.hartwig.hmftools.knowledgebaseimporter.Knowledgebase
import com.hartwig.hmftools.knowledgebaseimporter.civic.input.CivicEvidenceInput
import com.hartwig.hmftools.knowledgebaseimporter.civic.input.CivicVariantInput
import com.hartwig.hmftools.knowledgebaseimporter.diseaseOntology.DiseaseOntology
import com.hartwig.hmftools.knowledgebaseimporter.knowledgebases.ActionableRecord
import com.hartwig.hmftools.knowledgebaseimporter.knowledgebases.KnowledgebaseSource
import com.hartwig.hmftools.knowledgebaseimporter.knowledgebases.RecordAnalyzer
import com.hartwig.hmftools.knowledgebaseimporter.output.*

class Civic(variantsLocation: String, evidenceLocation: String, diseaseOntology: DiseaseOntology,
            private val recordAnalyzer: RecordAnalyzer, treatmentTypeMap: Map<String, String>) :
        Knowledgebase, KnowledgebaseSource<CivicRecord, ActionableRecord> {

    // KODU: This is a TP53 B-level evidence item that is dubious, so is filtered out.
    // KODU: See also https://civicdb.org/events/genes/45/summary/variants/222/summary/evidence/1481/summary#evidence
    private val blacklistedEvidenceIds = setOf("1481")

    override val source = "civic"
    override val knownVariants by lazy { recordAnalyzer.knownVariants(listOf(this)).distinct() }
    override val knownFusionPairs by lazy { knownKbRecords.flatMap { it.events }.filterIsInstance<FusionPair>().distinct() }
    override val promiscuousGenes by lazy { knownKbRecords.flatMap { it.events }.filterIsInstance<PromiscuousGene>().distinct() }
    override val actionableVariants by lazy { actionableKbItems.filterIsInstance<ActionableVariantOutput>() }
    override val actionableCNVs by lazy { actionableKbItems.filterIsInstance<ActionableCNVOutput>() }
    override val actionableFusionPairs by lazy { actionableKbItems.filterIsInstance<ActionableFusionPairOutput>() }
    override val actionablePromiscuousGenes by lazy { actionableKbItems.filterIsInstance<ActionablePromiscuousGeneOutput>() }
    override val actionableRanges by lazy { actionableKbItems.filterIsInstance<ActionableGenomicRangeOutput>() }
    override val cancerTypes by lazy {
        actionableKbRecords.flatMap { it.cancerDoids.entries }
                .associateBy({ it.key }, { diseaseOntology.findDoids(it.key) + diseaseOntology.findDoids(it.value) })
    }

    override val knownKbRecords by lazy { preProcessCivicRecords(variantsLocation, evidenceLocation, treatmentTypeMap) }
    override val actionableKbRecords by lazy { knownKbRecords }
    private val actionableKbItems by lazy { recordAnalyzer.actionableItems(listOf(this)).distinct() }

    private fun preProcessCivicRecords(variantFileLocation: String, evidenceFileLocation: String,
                                       treatmentTypeMap: Map<String, String>): List<CivicRecord> {
        val variantEvidenceMap = readEvidenceMap(evidenceFileLocation, treatmentTypeMap)
        return CsvReader.readTSVByName<CivicVariantInput>(variantFileLocation).mapNotNull { it.corrected() }.map {
            CivicRecord(it, variantEvidenceMap[it.variant_id])
        }
    }

    private fun readEvidenceMap(evidenceLocation: String, treatmentTypeMap: Map<String, String>): Multimap<String, CivicEvidence> {
        val civicApi = CivicApiWrapper()
        val drugInteractionMap = civicApi.drugInteractionMap
        val evidenceMap = ArrayListMultimap.create<String, CivicEvidence>()
        CsvReader.readTSVByName<CivicEvidenceInput>(evidenceLocation).mapNotNull { it.corrected() }.map {
            if (!blacklistedEvidenceIds.contains(it.evidence_id)) {
                evidenceMap.put(it.variant_id, CivicEvidence(it, drugInteractionMap, treatmentTypeMap))
            }
        }
        civicApi.releaseResources()
        return evidenceMap
    }
}
