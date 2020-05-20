package com.hartwig.hmftools.gripss.store

import com.hartwig.hmftools.gripss.StructuralVariantContext
import java.util.*
import kotlin.collections.HashMap

class VariantStore(
        private val variants: List<StructuralVariantContext>,
        private val variantsById: Map<String, StructuralVariantContext>,
        private val variantsByChromosome: Map<String, List<StructuralVariantContext>>) {

    companion object Factory {
        operator fun invoke(variants: List<StructuralVariantContext>): VariantStore {
            val variantsById = HashMap<String, StructuralVariantContext>()
            val variantsByChromosome = HashMap<String, MutableList<StructuralVariantContext>>()
            for (variant in variants) {
                variantsById[variant.vcfId] = variant
                variantsByChromosome.computeIfAbsent(variant.contig) { mutableListOf() }.add(variant)
            }

            return VariantStore(variants, variantsById, variantsByChromosome)
        }
    }

    fun selectAll(): List<StructuralVariantContext> {
        return variants
    }

    fun select(vcfId: String): StructuralVariantContext {
        return variantsById[vcfId]!!
    }

    fun contigsWithVariants(): Set<String> {
        return variantsByChromosome.keys
    }

    fun selectContigVariants(contig: String): List<StructuralVariantContext> {
        return variantsByChromosome.getOrDefault(contig, Collections.emptyList())
    }

    fun selectNearbyJoined(variant: StructuralVariantContext, maxDistance: Int): List<StructuralVariantContext> {
        val leftFilter = {other: StructuralVariantContext -> other.minStart <= variant.maxStart}
        val rightFilter = {other: StructuralVariantContext -> other.maxStart >= variant.minStart}
        val directionFilter = if (variant.orientation == 1.toByte()) leftFilter else rightFilter

        return selectNearby(variant, maxDistance, directionFilter)
    }

    fun selectNearby(variant: StructuralVariantContext, maxDistance: Int, filter: (StructuralVariantContext) -> Boolean = { _ -> true }): List<StructuralVariantContext> {
        val minStart = variant.minStart - maxDistance
        val maxStart = variant.maxStart + maxDistance
        val idFilter = {other: StructuralVariantContext -> variant.vcfId != other.vcfId && variant.mateId?.equals(other.vcfId) != true}
        val overlapFilter = { other: StructuralVariantContext -> other.minStart <= maxStart && other.maxStart >= minStart }
        val orientationFilter = { other: StructuralVariantContext -> other.orientation != variant.orientation }

        return variantsByChromosome.getOrDefault(variant.contig, Collections.emptyList()).filter { x -> idFilter(x) && overlapFilter(x) && orientationFilter(x) && filter(x) }
    }


    fun selectAlternatives(variant: StructuralVariantContext): List<StructuralVariantContext> {
        return selectOthersInConfidenceIntervals(variant) { x -> x.mateId?.equals(variant.vcfId) != true && x.orientation == variant.orientation }
    }

    fun selectOthersInConfidenceIntervals(variant: StructuralVariantContext, filter: (StructuralVariantContext) -> Boolean = { true }): List<StructuralVariantContext> {
        return variantsByChromosome.getOrDefault(variant.contig, Collections.emptyList()).filter { other ->
            return@filter filter(other) && variant.confidenceIntervalsOverlap(other) && other.vcfId != variant.vcfId
        }
    }


}