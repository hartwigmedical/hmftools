package com.hartwig.hmftools.knowledgebaseimporter

import com.hartwig.hmftools.knowledgebaseimporter.output.ActionableCNVOutput
import com.hartwig.hmftools.knowledgebaseimporter.output.ActionableFusionOutput
import com.hartwig.hmftools.knowledgebaseimporter.output.ActionableVariantOutput
import com.hartwig.hmftools.knowledgebaseimporter.output.KnownVariantOutput

interface Knowledgebase {
    val knownVariants: List<KnownVariantOutput>
    val knownFusionPairs: List<FusionPair>
    val promiscuousGenes: List<PromiscuousGene>
    val actionableVariants: List<ActionableVariantOutput>
    val actionableCNVs: List<ActionableCNVOutput>
    val actionableFusions: List<ActionableFusionOutput>
}
