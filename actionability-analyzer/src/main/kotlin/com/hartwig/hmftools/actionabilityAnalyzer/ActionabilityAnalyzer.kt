package com.hartwig.hmftools.actionabilityAnalyzer

import com.hartwig.hmftools.common.copynumber.CopyNumberAlteration
import com.hartwig.hmftools.knowledgebaseimporter.output.*
import com.hartwig.hmftools.knowledgebaseimporter.readCSVRecords
import com.hartwig.hmftools.knowledgebaseimporter.readTSVRecords
import com.hartwig.hmftools.patientdb.data.PotentialActionableCNV
import com.hartwig.hmftools.patientdb.data.PotentialActionableFusion
import com.hartwig.hmftools.patientdb.data.PotentialActionableVariant
import org.apache.commons.csv.CSVRecord

class ActionabilityAnalyzer(actionableVariantsLocation: String, fusionPairsLocation: String, promiscuousFiveLocation: String,
                            promiscuousThreeLocation: String, cnvsLocation: String, cancerTypeLocation: String) {
    companion object {
        private fun <T> createActionabilityMap(items: List<ActionableItem<T>>): Map<T, List<ActionableTreatment>> {
            return items.groupBy { it.event }.mapValues { (_, actionableOutputs) ->
                ActionableTreatment(actionableOutputs)
            }
        }

        private fun readActionableVariants(fileLocation: String): List<ActionableVariantOutput> {
            return readTSVRecords(fileLocation) {
                ActionableVariantOutput(it["gene"], SomaticVariantEvent(it["gene"], it["chromosome"], it["position"], it["ref"], it["alt"]),
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
                ActionableCNVOutput(CnvEvent(it["gene"], it["cnvType"]), readActionability(it))
            }
        }

        private fun readActionability(record: CSVRecord): Actionability {
            return Actionability(record["source"], record["cancerType"], record["drug"].orEmpty(), record["level"],
                                 record["significance"].orEmpty(), record["evidenceType"].orEmpty(), HmfLevel.valueOf(record["hmfLevel"]),
                                 HmfResponse.valueOf(record["hmfResponse"]))
        }

        private fun readCancerTypeMapping(fileLocation: String): Map<String, Set<String>> {
            return readTSVRecords(fileLocation) { CancerTypeDoidOutput(it["cancerType"], it["doids"].orEmpty()) }
                    .map { Pair(it.cancerType, it.doidSet.split(";").filterNot { it.isBlank() }.toSet()) }
                    .toMap()
        }

        private fun readPrimaryTumorMapping(): Map<String, Set<String>> {
            return readCSVRecords(this::class.java.getResourceAsStream("/primary_tumor_locations.csv")) {
                Pair(it["primaryTumorLocation"], it["doids"].orEmpty().split(";").filterNot { it.isBlank() }.toSet())
            }.toMap()
        }
    }

    private val variantActionabilityMap = createActionabilityMap(readActionableVariants(actionableVariantsLocation))
    private val fusionActionabilityMap = createActionabilityMap(readActionableFusionPairs(fusionPairsLocation))
    private val promiscuousFiveActionabilityMap = createActionabilityMap(readActionablePromiscuousGenes(promiscuousFiveLocation))
    private val promiscuousThreeActionabilityMap = createActionabilityMap(readActionablePromiscuousGenes(promiscuousThreeLocation))
    private val cnvActionabilityMap = createActionabilityMap(readActionableCNVs(cnvsLocation))
    private val cancerTypeMapping = readCancerTypeMapping(cancerTypeLocation)
    private val tumorLocationMapping = readPrimaryTumorMapping()

    fun actionabilityForVariant(variant: PotentialActionableVariant): Set<ActionabilityOutput> {
        val variantKey = SomaticVariantEvent(variant.gene(), variant.chromosome(), variant.position().toString(), variant.ref(),
                                             variant.alt())
        return getActionability(variantActionabilityMap, variantKey, variant.sampleId(), variant.primaryTumorLocation()).toSet()
    }

    fun actionabilityForFusion(fusion: PotentialActionableFusion): Set<ActionabilityOutput> {
        val fiveGene = fusion.fiveGene()
        val threeGene = fusion.threeGene()
        val sampleId = fusion.sampleId()
        val cancerType = fusion.primaryTumorLocation()
        val fusionPairActionability = getActionability(fusionActionabilityMap, FusionPair(fiveGene, threeGene), sampleId, cancerType)
        val promiscuousFiveActionability = getActionability(promiscuousFiveActionabilityMap, PromiscuousGene(fiveGene), sampleId,
                                                            cancerType)
        val promiscuousThreeActionability = getActionability(promiscuousThreeActionabilityMap, PromiscuousGene(threeGene), sampleId,
                                                             cancerType)
        return (fusionPairActionability + promiscuousFiveActionability + promiscuousThreeActionability).toSet()
    }

    fun actionabilityForCNV(cnv: PotentialActionableCNV): Set<ActionabilityOutput> {
        val cnvType = if (cnv.alteration() == CopyNumberAlteration.GAIN) "Amplification" else "Deletion"
        val cnvEvent = CnvEvent(cnv.gene(), cnvType)
        return getActionability(cnvActionabilityMap, cnvEvent, cnv.sampleId(), cnv.primaryTumorLocation()).toSet()
    }

    private fun <T> getActionability(actionabilityMap: Map<T, List<ActionableTreatment>>, event: T,
                                     sampleId: String, cancerType: String?): List<ActionabilityOutput> {
        return actionabilityMap[event].orEmpty().map {
            val treatmentType = getTreatmentType(cancerType, it.actionability.cancerType)
            ActionabilityOutput(sampleId, cancerType, treatmentType, it)
        }
    }

    private fun getTreatmentType(patientCancerType: String?, treatmentCancerType: String): String {
        val cancerDoids = tumorLocationMapping[patientCancerType].orEmpty()
        if (cancerDoids.isEmpty()) return "Unknown"
        val diseaseDoids = cancerTypeMapping[treatmentCancerType].orEmpty()
        if (diseaseDoids.isEmpty()) return "Unknown"
        val match = cancerDoids.any { diseaseDoids.contains(it) }
        return if (match) "On-label" else "Off-label"
    }
}
