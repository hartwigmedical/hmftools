package com.hartwig.hmftools.knowledgebaseimporter

import com.hartwig.hmftools.knowledgebaseimporter.output.*

interface Knowledgebase {
    val source: String
    val knownVariants: List<KnownVariantOutput>
    val knownFusionPairs: List<FusionPair>
    val promiscuousGenes: List<PromiscuousGene>
    val actionableVariants: List<ActionableVariantOutput>
    val actionableCNVs: List<ActionableCNVOutput>
    val actionableFusions: List<ActionableFusionOutput>
}
