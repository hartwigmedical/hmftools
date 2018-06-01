package com.hartwig.hmftools.knowledgebaseimporter

import com.hartwig.hmftools.knowledgebaseimporter.output.*

interface Knowledgebase {
    val source: String
    val knownVariants: List<KnownVariantOutput>
    val knownFusionPairs: List<FusionPair>
    val promiscuousGenes: List<PromiscuousGene>
    val actionableVariants: List<ActionableVariantOutput>
    val actionableCNVs: List<ActionableCNVOutput>
    val actionableFusionPairs: List<ActionableFusionPairOutput>
    val actionablePromiscuousGenes: List<ActionablePromiscuousGeneOutput>
    val cancerTypes: Map<String, Set<String>>
}
