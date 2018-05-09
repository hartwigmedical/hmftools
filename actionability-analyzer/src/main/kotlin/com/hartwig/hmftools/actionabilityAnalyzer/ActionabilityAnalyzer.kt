package com.hartwig.hmftools.actionabilityAnalyzer

import com.hartwig.hmftools.knowledgebaseimporter.output.*
import com.hartwig.hmftools.knowledgebaseimporter.readTSVRecords
import com.hartwig.hmftools.patientdb.data.PotentialActionableVariant
import org.apache.commons.csv.CSVRecord

class ActionabilityAnalyzer(actionableVariantsLocation: String, fustionPairsLocation: String, promiscuousFiveLocation: String,
                            promiscuousThreeLocation: String, cnvsLocation: String) {
    companion object {
        private fun <T> createActionabilityMap(items: List<ActionableItem<T>>): Map<T, List<ActionableTreatment>> {
            return items.groupBy { it.event }.mapValues { (_, actionableOutputs) ->
                ActionableTreatment(actionableOutputs)
            }
        }

        private fun readActionableVariants(fileLocation: String): List<ActionableVariantOutput> {
            return readTSVRecords(fileLocation) {
                ActionableVariantOutput(it["gene"], SomaticVariantEvent(it["chromosome"], it["position"], it["ref"], it["alt"]),
                                        readActionability(it))
            }
        }

        private fun readActionableFusionPairs(fileLocation: String): List<ActionableFusionOutput> {
            return readTSVRecords(fileLocation) {
                ActionableFusionOutput(FusionPair(it["fiveGene"], it["threeGene"]), readActionability(it))
            }
        }

        private fun readActionablePromiscuousGenes(fileLocation: String): List<ActionableFusionOutput> {
            return readTSVRecords(fileLocation) {
                ActionableFusionOutput(PromiscuousGene(it["gene"]), readActionability(it))
            }
        }

        private fun readActionableCNVs(fileLocation: String): List<ActionableCNVOutput> {
            return readTSVRecords(fileLocation) {
                ActionableCNVOutput(CnvEvent(it["gene"], it["type"]), readActionability(it))
            }
        }

        private fun readActionability(record: CSVRecord): Actionability {
            return Actionability(record["source"], record["cancerType"], record["drug"].orEmpty(), record["level"],
                                 record["significance"].orEmpty(), record["evidenceType"].orEmpty())
        }

        private fun <T> getActionability(actionabilityMap: Map<T, List<ActionableTreatment>>, event: T,
                                         sampleId: String): List<ActionabilityOutput> {
            return actionabilityMap[event].orEmpty().map { ActionabilityOutput(sampleId, it) }
        }
    }

    private val variantActionabilityMap = createActionabilityMap(readActionableVariants(actionableVariantsLocation))
    private val fusionActionabilityMap = createActionabilityMap(readActionableFusionPairs(fustionPairsLocation))
    private val promiscuousFiveActionabilityMap = createActionabilityMap(readActionablePromiscuousGenes(promiscuousFiveLocation))
    private val promiscuousThreeActionabilityMap = createActionabilityMap(readActionablePromiscuousGenes(promiscuousThreeLocation))
    private val cnvActionabilityMap = createActionabilityMap(readActionableCNVs(cnvsLocation))

    fun actionabilityForVariant(variant: PotentialActionableVariant): Set<ActionabilityOutput> {
        val variantKey = SomaticVariantEvent(variant.chromosome(), variant.position().toString(), variant.ref(), variant.alt())
        return getActionability(variantActionabilityMap, variantKey, variant.sampleId()).toSet()
    }

    fun actionabilityForFusion(sampleId: String, fiveGene: String, threeGene: String): Set<ActionabilityOutput> {
        val fusionPairActionability = getActionability(fusionActionabilityMap, FusionPair(fiveGene, threeGene), sampleId)
        val promiscuousFiveActionability = getActionability(promiscuousFiveActionabilityMap, PromiscuousGene(fiveGene), sampleId)
        val promiscuousThreeActionability = getActionability(promiscuousThreeActionabilityMap, PromiscuousGene(threeGene), sampleId)
        return (fusionPairActionability + promiscuousFiveActionability + promiscuousThreeActionability).toSet()
    }

    fun actionabilityForCNV(sampleId: String, gene: String, cnvType: String): Set<ActionabilityOutput> {
        val cnvEvent = CnvEvent(gene, cnvType)
        return getActionability(cnvActionabilityMap, cnvEvent, sampleId).toSet()
    }
}
