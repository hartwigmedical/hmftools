package com.hartwig.hmftools.knowledgebaseimporter.output

import com.hartwig.hmftools.extensions.csv.CsvData
import com.hartwig.hmftools.knowledgebaseimporter.knowledgebases.ActionableItem

sealed class CsvOutput<T : ActionableEvent> : ActionableItem<T>, CsvData

data class ActionableCNVOutput(override val event: CnvEvent, override val actionability: Actionability) : CsvOutput<CnvEvent>()

data class ActionableFusionPairOutput(override val event: FusionPair, override val actionability: Actionability) : CsvOutput<FusionPair>()

data class ActionablePromiscuousGeneOutput(override val event: PromiscuousGene,
                                           override val actionability: Actionability) : CsvOutput<PromiscuousGene>()

data class ActionableVariantOutput(override val event: SomaticVariantEvent,
                                   override val actionability: Actionability) : CsvOutput<SomaticVariantEvent>()
