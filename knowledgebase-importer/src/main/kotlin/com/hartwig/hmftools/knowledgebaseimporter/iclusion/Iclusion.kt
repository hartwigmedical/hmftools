package com.hartwig.hmftools.knowledgebaseimporter.iclusion

import com.hartwig.hmftools.apiclients.iclusion.data.IclusionStudyDetails
import com.hartwig.hmftools.knowledgebaseimporter.Knowledgebase
import com.hartwig.hmftools.knowledgebaseimporter.diseaseOntology.DiseaseOntology
import com.hartwig.hmftools.knowledgebaseimporter.gene.GeneModel
import com.hartwig.hmftools.knowledgebaseimporter.knowledgebases.ActionableRecord
import com.hartwig.hmftools.knowledgebaseimporter.knowledgebases.KnowledgebaseSource
import com.hartwig.hmftools.knowledgebaseimporter.knowledgebases.RecordAnalyzer
import com.hartwig.hmftools.knowledgebaseimporter.output.*

class Iclusion(iclusionStudies: List<IclusionStudyDetails>, diseaseOntology: DiseaseOntology, recordAnalyzer: RecordAnalyzer) :
        Knowledgebase, KnowledgebaseSource<IclusionRecord, ActionableRecord> {
    override val source: String = "iclusion"

    private val geneModel = GeneModel()
    private val geneToTranscriptMap = iclusionStudies.flatMap { it -> it.mutations.map { it.geneName } }.distinct().map {
        Pair(it, geneModel.hmfTranscriptForGene(it))
    }.toMap()

    override val knownVariants: List<KnownVariantOutput> = listOf()
    override val knownFusionPairs: List<FusionPair> = listOf()
    override val promiscuousGenes: List<PromiscuousGene> = listOf()

    override val actionableVariants by lazy { actionableKbItems.filterIsInstance<ActionableVariantOutput>() }
    override val actionableCNVs by lazy { actionableKbItems.filterIsInstance<ActionableCNVOutput>() }
    override val actionableFusionPairs by lazy { actionableKbItems.filterIsInstance<ActionableFusionPairOutput>() }
    override val actionablePromiscuousGenes by lazy { actionableKbItems.filterIsInstance<ActionablePromiscuousGeneOutput>() }
    override val actionableRanges by lazy { actionableKbItems.filterIsInstance<ActionableGenomicRangeOutput>() }
    override val cancerTypes by lazy {
        actionableKbRecords.flatMap { it.doids.entries }
                .associateBy({ it.key }, { it -> it.value.flatMap { diseaseOntology.findDoids(it) }.toSet() })
    }

    override val knownKbRecords: List<IclusionRecord> by lazy { iclusionStudies.flatMap { IclusionRecord(it, geneToTranscriptMap) } }

    override val actionableKbRecords by lazy { knownKbRecords }
    private val actionableKbItems by lazy { recordAnalyzer.actionableItems(listOf(this)).distinct() }
}
