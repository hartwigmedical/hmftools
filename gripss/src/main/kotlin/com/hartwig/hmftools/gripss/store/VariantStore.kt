package com.hartwig.hmftools.gripss.store

import com.hartwig.hmftools.gripss.StructuralVariantContext
import java.util.*
import kotlin.collections.HashMap
import kotlin.math.abs
import kotlin.math.max

class VariantStore(
        private val variants: List<StructuralVariantContext>,
        private val variantIndexesById: Map<String, Int>,
        private val variantsByChromosome: Map<String, List<StructuralVariantContext>>) {

    companion object Factory {
        operator fun invoke(variants: List<StructuralVariantContext>): VariantStore {
            val variantIndexesById = HashMap<String, Int>()
            val variantsByChromosome = HashMap<String, MutableList<StructuralVariantContext>>()
            for (i in variants.indices) {
                val variant = variants[i]
                variantIndexesById[variant.vcfId] = i
                variantsByChromosome.computeIfAbsent(variant.contig) { mutableListOf() }.add(variant)
            }

            return VariantStore(variants, variantIndexesById, variantsByChromosome)
        }
    }

    fun selectAll(): List<StructuralVariantContext> {
        return variants
    }

    fun selectAllInContig(contig: String): List<StructuralVariantContext> {
        return variantsByChromosome.getOrDefault(contig, Collections.emptyList())
    }

    fun select(vcfId: String): StructuralVariantContext {
        return variants[variantIndexesById[vcfId]!!]
    }

    fun selectOthersNearby(variant: StructuralVariantContext, maxDistance: Pair<Int, Int>, filter: (StructuralVariantContext) -> Boolean = { _ -> true }): List<StructuralVariantContext> {
        val minStart = variant.minStart - maxDistance.first
        val maxStart = variant.maxStart + maxDistance.second
        val idFilter = { other: StructuralVariantContext -> variant.vcfId != other.vcfId && variant.mateId?.equals(other.vcfId) != true }
        val overlapFilter = { other: StructuralVariantContext -> other.minStart <= maxStart && other.maxStart >= minStart }
        return variantsByChromosome.getOrDefault(variant.contig, Collections.emptyList()).filter { x -> idFilter(x) && overlapFilter(x) && filter(x) }
    }

    @Deprecated("Use other", ReplaceWith("selectOthersNearby(variant, Pair(maxDistance, maxDistance), filter)"))
    fun selectOthersNearby(variant: StructuralVariantContext, maxDistance: Int, filter: (StructuralVariantContext) -> Boolean = { _ -> true }): List<StructuralVariantContext> {
        return selectOthersNearby(variant, Pair(maxDistance, maxDistance), filter)
    }

    fun selectOthersNearby(variant: StructuralVariantContext, additionalDistance: Pair<Int, Int>, seekDistance: Int, filter: (StructuralVariantContext) -> Boolean): List<StructuralVariantContext> {
        return selectOtherIndicesNearby(variantIndexesById[variant.vcfId]!!, seekDistance, additionalDistance, filter)
    }

    private fun selectOtherIndicesNearby(i: Int, maxSeekDistance: Int, maxAdditionalDistance: Pair<Int, Int>, filter: (StructuralVariantContext) -> Boolean): List<StructuralVariantContext> {
        val variant = variants[i]
        val minStart = variant.minStart - abs(maxAdditionalDistance.first)
        val maxStart = variant.maxStart + maxAdditionalDistance.second
        val idFilter = { other: StructuralVariantContext -> variant.vcfId != other.vcfId && variant.mateId?.equals(other.vcfId) != true }
        val overlapFilter = { other: StructuralVariantContext -> other.minStart <= maxStart && other.maxStart >= minStart }
        val allFilters = { other: StructuralVariantContext -> idFilter(other) && overlapFilter(other) && filter(other) }

        val result = mutableListOf<StructuralVariantContext>()
        // Look forwards
        for (j in i + 1 until variants.size) {
            val other = variants[j]
            if (other.minStart > variant.maxStart + maxSeekDistance) {
                break
            } else if (allFilters(other)) {
                result.add(other)
            }
        }

        // Look backwards
        for (j in max(0, i - 1) downTo 0) {
            val other = variants[j]
            if (other.maxStart < variant.minStart - maxSeekDistance) {
                break
            } else if (allFilters(other)) {
                result.add(other)
            }
        }

        return result
    }

}