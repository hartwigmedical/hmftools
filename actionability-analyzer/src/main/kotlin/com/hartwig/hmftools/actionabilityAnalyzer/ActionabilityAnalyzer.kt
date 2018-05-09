package com.hartwig.hmftools.actionabilityAnalyzer

import com.hartwig.hmftools.knowledgebaseimporter.output.Actionability
import com.hartwig.hmftools.knowledgebaseimporter.output.ActionableItem
import com.hartwig.hmftools.knowledgebaseimporter.output.ActionableVariantOutput
import com.hartwig.hmftools.knowledgebaseimporter.output.SomaticVariantEvent
import com.hartwig.hmftools.knowledgebaseimporter.readTSVRecords
import com.hartwig.hmftools.patientdb.data.PotentialActionableVariant

class ActionabilityAnalyzer(actionableVariantsLocation: String) {
    companion object {
        private fun <T> createActionabilityMap(items: List<ActionableItem<T>>): Map<T, List<ActionableTreatment>> {
            return items.groupBy { it.event }.mapValues { (_, actionableOutputs) ->
                ActionableTreatment(actionableOutputs)
            }
        }

        private fun readActionableVariants(actionableVariantsLocation: String): List<ActionableVariantOutput> {
            return readTSVRecords(actionableVariantsLocation) {
                ActionableVariantOutput(it["gene"], SomaticVariantEvent(it["chromosome"], it["position"], it["ref"], it["alt"]),
                                        Actionability(it["source"], it["cancerType"], it["drug"] ?: "", it["level"],
                                                      it["significance"] ?: "", it["evidenceType"] ?: ""))
            }
        }
    }

    private val variantActionabilityMap = createActionabilityMap(readActionableVariants(actionableVariantsLocation))
    //    private val fusionActionabilityMap: String = TODO("not implemented")
    //    private val cnvActionabilityMap: String = TODO("not implemented")

    fun actionabilityForVariant(variant: PotentialActionableVariant): List<ActionabilityOutput> {
        val variantKey = SomaticVariantEvent(variant.chromosome(), variant.position().toString(), variant.ref(), variant.alt())
        return variantActionabilityMap[variantKey]?.map { ActionabilityOutput(variant.sampleId(), it) } ?: listOf()
    }
}
