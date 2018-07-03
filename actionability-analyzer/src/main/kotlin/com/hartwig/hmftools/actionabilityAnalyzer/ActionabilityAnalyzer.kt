package com.hartwig.hmftools.actionabilityAnalyzer

import com.github.davidmoten.rtree.RTree
import com.github.davidmoten.rtree.geometry.Geometries
import com.github.davidmoten.rtree.geometry.Line
import com.google.common.annotations.VisibleForTesting
import com.hartwig.hmftools.common.copynumber.CopyNumberAlteration
import com.hartwig.hmftools.extensions.csv.CsvReader
import com.hartwig.hmftools.knowledgebaseimporter.knowledgebases.ActionableItem
import com.hartwig.hmftools.knowledgebaseimporter.output.*
import com.hartwig.hmftools.knowledgebaseimporter.readCSVRecords
import com.hartwig.hmftools.knowledgebaseimporter.readTSVRecords
import com.hartwig.hmftools.patientdb.data.PotentialActionableCNV
import com.hartwig.hmftools.patientdb.data.PotentialActionableFusion

class ActionabilityAnalyzer(private val sampleTumorLocationMap: Map<String, String>, actionableVariantsLocation: String,
                            fusionPairsLocation: String, promiscuousFiveLocation: String, promiscuousThreeLocation: String,
                            cnvsLocation: String, cancerTypeLocation: String, actionableRangesLocation: String) {
    companion object {
        private fun <T : ActionableEvent, R> createActionabilityMap(items: List<ActionableItem<T>>,
                                                                    keyMapper: (T) -> R): Map<R, List<ActionableTreatment>> {
            return items.groupBy { keyMapper(it.event) }.mapValues { (_, actionableOutputs) -> ActionableTreatment(actionableOutputs) }
        }

        private fun createActionabilityTree(
                items: List<ActionableGenomicRangeOutput>): Map<String, RTree<List<ActionableTreatment>, Line>> {
            return items.groupBy { it.event.chromosome }.mapValues { rangesToRtree(it.value) }
        }

        private fun rangesToRtree(items: List<ActionableGenomicRangeOutput>): RTree<List<ActionableTreatment>, Line> {
            return items.groupBy { it.event }.entries.fold(RTree.create<List<ActionableTreatment>, Line>()) { rtree, entry ->
                val treatments = ActionableTreatment(entry.value)
                val range = entry.key
                rtree.add(treatments, genomicPositionsToLine(range.start.toInt(), range.stop.toInt()))
            }
        }

        private fun genomicPositionsToLine(start: Int, end: Int): Line {
            return Geometries.line(start.toDouble(), 0.toDouble(), end.toDouble(), 0.toDouble())
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
    private val rangeActionabilityTree = createActionabilityTree(CsvReader.readTSV(actionableRangesLocation))
    private val cancerTypeMapping = readCancerTypeMapping(cancerTypeLocation)
    private val tumorLocationMapping = primaryTumorMapping()

    fun actionabilityForVariant(variant: CohortMutation): Set<ActionabilityOutput> {
        val variantKey = VariantKey(variant)
        val cancerType = sampleTumorLocationMap[variant.sampleId]
        val variantString = potentialVariantString(variant)
        return getActionability(variantActionabilityMap, variantKey, variantString, variant.sampleId, cancerType, variant.gene, "",
                                variant.type).toSet()
    }

    fun rangeActionabilityForVariant(variant: CohortMutation): Set<ActionabilityOutput> {
        if (!variant.potentiallyActionable) return emptySet()
        val cancerType = sampleTumorLocationMap[variant.sampleId]
        val variantPosition = genomicPositionsToLine(variant.position.toInt(), variant.position.toInt() + variant.ref.length - 1)
        val searchResults = rangeActionabilityTree[variant.chromosome]?.search(variantPosition.mbr(), 0.1)
        searchResults ?: return emptySet()
        val treatments = searchResults.toBlocking().toIterable().flatMap { it.value() }
        val variantString = potentialVariantString(variant)
        return treatments.map {
            val treatmentType = getTreatmentType(cancerType, it.actionability.cancerType)
            ActionabilityOutput(variant.sampleId, variantString, variant.gene, "", variant.type, cancerType, treatmentType, it)
        }.toSet()
    }

    fun actionabilityForFusion(fusion: PotentialActionableFusion): Set<ActionabilityOutput> {
        val fiveGene = fusion.fiveGene()
        val threeGene = fusion.threeGene()
        val sampleId = fusion.sampleId()
        val cancerType = sampleTumorLocationMap[fusion.sampleId()]
        val eventString = "$fiveGene - $threeGene fusion"
        val fusionPairActionability = getActionability(fusionActionabilityMap, FusionPair(fiveGene, threeGene), eventString, sampleId,
                                                       cancerType, fusion.fiveGene(), fusion.threeGene(), "Fusion")
        val promiscuousFiveActionability = getActionability(promiscuousFiveActionabilityMap, PromiscuousGene(fiveGene), eventString,
                                                            sampleId, cancerType, fusion.fiveGene(), fusion.threeGene(), "Fusion")
        val promiscuousThreeActionability = getActionability(promiscuousThreeActionabilityMap, PromiscuousGene(threeGene), eventString,
                                                             sampleId, cancerType, fusion.fiveGene(), fusion.threeGene(), "Fusion")
        return (fusionPairActionability + promiscuousFiveActionability + promiscuousThreeActionability).toSet()
    }

    fun actionabilityForCNV(cnv: PotentialActionableCNV): Set<ActionabilityOutput> {
        val cnvType = if (cnv.alteration() == CopyNumberAlteration.GAIN) "Amplification" else "Deletion"
        val cnvEvent = CnvEvent(cnv.gene(), cnvType)
        return getActionability(cnvActionabilityMap, cnvEvent, cnvEvent.eventString(), cnv.sampleId(),
                                sampleTumorLocationMap[cnv.sampleId()], cnv.gene(), "", cnvType).toSet()
    }

    private fun <T> getActionability(actionabilityMap: Map<T, List<ActionableTreatment>>, event: T, eventString: String,
                                     sampleId: String, cancerType: String?, gene: String, partner: String,
                                     eventType: String): List<ActionabilityOutput> {
        return actionabilityMap[event].orEmpty().map {
            val treatmentType = getTreatmentType(cancerType, it.actionability.cancerType)
            ActionabilityOutput(sampleId, eventString, gene, partner, eventType, cancerType, treatmentType, it)
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

    private fun potentialVariantString(variant: CohortMutation): String {
        return "${variant.impact} ${variant.chromosome} ${variant.position}: ${variant.ref} -> ${variant.alt}"
    }
}
