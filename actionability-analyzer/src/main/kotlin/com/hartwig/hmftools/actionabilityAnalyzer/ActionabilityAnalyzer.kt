package com.hartwig.hmftools.actionabilityAnalyzer

import com.google.common.annotations.VisibleForTesting
import com.hartwig.hmftools.common.copynumber.CopyNumberAlteration
import com.hartwig.hmftools.extensions.csv.CsvReader
import com.hartwig.hmftools.knowledgebaseimporter.knowledgebases.ActionableItem
import com.hartwig.hmftools.knowledgebaseimporter.output.*
import com.hartwig.hmftools.knowledgebaseimporter.readCSVRecords
import com.hartwig.hmftools.knowledgebaseimporter.readTSVRecords
import com.hartwig.hmftools.patientdb.data.PotentialActionableCNV
import com.hartwig.hmftools.patientdb.data.PotentialActionableFusion
import com.hartwig.hmftools.patientdb.data.PotentialActionableVariant

class ActionabilityAnalyzer(private val sampleTumorLocationMap: Map<String, String>, actionableVariantsLocation: String,
                            fusionPairsLocation: String, promiscuousFiveLocation: String, promiscuousThreeLocation: String,
                            cnvsLocation: String, cancerTypeLocation: String) {
    companion object {
        private fun <T : ActionableEvent, R> createActionabilityMap(items: List<ActionableItem<T>>,
                                                                    keyMapper: (T) -> R): Map<R, List<ActionableTreatment>> {
            return items.groupBy { keyMapper(it.event) }.mapValues { (_, actionableOutputs) -> ActionableTreatment(actionableOutputs) }
        }

        private fun readCancerTypeMapping(fileLocation: String): Map<String, Set<String>> {
            return readTSVRecords(fileLocation) { CancerTypeDoidOutput(it["cancerType"], it["doids"].orEmpty()) }
                    .map { Pair(it.cancerType, it.doidSet.split(";").filterNot { it.isBlank() }.toSet()) }
                    .toMap()
        }

        @VisibleForTesting
        fun primaryTumorMapping(): Map<String, Set<String>> {
            return readCSVRecords(this::class.java.getResourceAsStream("/primary_tumor_locations.csv")) {
                Pair(it["primaryTumorLocation"], it["doids"].orEmpty().split(";").filterNot { it.isBlank() }.toSet())
            }.toMap()
        }
    }

    private val variantActionabilityMap =
            createActionabilityMap(CsvReader.readTSV<ActionableVariantOutput>(actionableVariantsLocation)) { VariantKey(it) }
    private val fusionActionabilityMap =
            createActionabilityMap(CsvReader.readTSV<ActionableFusionPairOutput>(fusionPairsLocation)) { it }
    private val promiscuousFiveActionabilityMap =
            createActionabilityMap(CsvReader.readTSV<ActionablePromiscuousGeneOutput>(promiscuousFiveLocation)) { it }
    private val promiscuousThreeActionabilityMap =
            createActionabilityMap(CsvReader.readTSV<ActionablePromiscuousGeneOutput>(promiscuousThreeLocation)) { it }
    private val cnvActionabilityMap = createActionabilityMap(CsvReader.readTSV<ActionableCNVOutput>(cnvsLocation)) { it }
    private val cancerTypeMapping = readCancerTypeMapping(cancerTypeLocation)
    private val tumorLocationMapping = primaryTumorMapping()

    fun actionabilityForVariant(variant: PotentialActionableVariant): Set<ActionabilityOutput> {
        val variantKey = VariantKey(variant)
        val cancerType = sampleTumorLocationMap[variant.sampleId()]
        return getActionability(variantActionabilityMap, variantKey, variant.sampleId(), cancerType).toSet()
    }

    fun actionabilityForFusion(fusion: PotentialActionableFusion): Set<ActionabilityOutput> {
        val fiveGene = fusion.fiveGene()
        val threeGene = fusion.threeGene()
        val sampleId = fusion.sampleId()
        val cancerType = sampleTumorLocationMap[fusion.sampleId()]
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
        return getActionability(cnvActionabilityMap, cnvEvent, cnv.sampleId(), sampleTumorLocationMap[cnv.sampleId()]).toSet()
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
