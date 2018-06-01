package com.hartwig.hmftools.knowledgebaseimporter.oncoKb

import com.hartwig.hmftools.knowledgebaseimporter.Knowledgebase
import com.hartwig.hmftools.knowledgebaseimporter.diseaseOntology.DiseaseOntology
import com.hartwig.hmftools.knowledgebaseimporter.knowledgebases.ActionableRecord
import com.hartwig.hmftools.knowledgebaseimporter.knowledgebases.KnowledgebaseSource
import com.hartwig.hmftools.knowledgebaseimporter.knowledgebases.RecordAnalyzer
import com.hartwig.hmftools.knowledgebaseimporter.output.*
import com.hartwig.hmftools.knowledgebaseimporter.readTSVRecords
import htsjdk.samtools.reference.IndexedFastaSequenceFile


class OncoKb(annotatedVariantsLocation: String, actionableVariantsLocation: String, transvarLocation: String,
             diseaseOntology: DiseaseOntology, private val reference: IndexedFastaSequenceFile,
             treatmentTypeMap: Map<String, String>) : Knowledgebase,
        KnowledgebaseSource<OncoKnownRecord, ActionableRecord> {
    override val source = "oncoKb"
    override val knownVariants by lazy { RecordAnalyzer(transvarLocation, reference).knownVariants(listOf(this)).distinct() }
    override val knownFusionPairs by lazy { knownKbRecords.flatMap { it.events }.filterIsInstance<FusionPair>().distinct() }
    override val promiscuousGenes by lazy { knownKbRecords.flatMap { it.events }.filterIsInstance<PromiscuousGene>().distinct() }
    override val actionableVariants by lazy { actionableKbItems.map { it.toActionableOutput() }.filterIsInstance<ActionableVariantOutput>() }
    override val actionableCNVs by lazy { actionableKbItems.map { it.toActionableOutput() }.filterIsInstance<ActionableCNVOutput>() }
    override val actionableFusionPairs by lazy { actionableKbItems.map { it.toActionableOutput() }.filterIsInstance<ActionableFusionPairOutput>() }
    override val actionablePromiscuousGenes by lazy { actionableKbItems.map { it.toActionableOutput() }.filterIsInstance<ActionablePromiscuousGeneOutput>() }
    override val cancerTypes by lazy {
        actionableKbRecords.flatMap { it.actionability }.map { it.cancerType }
                .associateBy({ it }, { diseaseOntology.findDoidsForCancerType(it) })
    }
    override val knownKbRecords by lazy { readTSVRecords(annotatedVariantsLocation) { OncoKnownRecord(it) } }
    override val actionableKbRecords by lazy { readTSVRecords(actionableVariantsLocation) { OncoActionableRecord(it, treatmentTypeMap) } }
    private val actionableKbItems by lazy { RecordAnalyzer(transvarLocation, reference).actionableItems(listOf(this)).distinct() }
}
