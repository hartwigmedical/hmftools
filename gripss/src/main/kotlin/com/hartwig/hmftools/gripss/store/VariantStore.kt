package com.hartwig.hmftools.gripss.store

import com.hartwig.hmftools.gripss.StructuralVariantContext
import java.util.*
import kotlin.collections.HashMap
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

    fun selectOthersNearby(variant: StructuralVariantContext, additionalDistance: Int, seekDistance: Int, filter: (StructuralVariantContext) -> Boolean): Collection<StructuralVariantContext> {
        return selectOthersNearby(variantIndexesById[variant.vcfId]!!, additionalDistance, seekDistance, filter)
    }

    private fun selectOthersNearby(i: Int, additionalDistance: Int,  maxSeekDistance: Int, filter: (StructuralVariantContext) -> Boolean): Collection<StructuralVariantContext> {
        val variant = variants[i]
        val minStart = variant.minStart - additionalDistance
        val maxStart = variant.maxStart + additionalDistance
        val idFilter = { other: StructuralVariantContext -> variant.vcfId != other.vcfId && variant.mateId?.equals(other.vcfId) != true }
        val overlapFilter = { other: StructuralVariantContext -> other.minStart <= maxStart && other.maxStart >= minStart }
        val allFilters = { other: StructuralVariantContext -> idFilter(other) && overlapFilter(other) && filter(other) }

        val result = ArrayDeque<StructuralVariantContext>()

        // Look backwards
        for (j in max(0, i - 1) downTo 0) {
            val other = variants[j]
            if (other.maxStart < variant.minStart - maxSeekDistance) {
                break
            } else if (allFilters(other)) {
                result.addFirst(other)
            }
        }

        // Look forwards
        for (j in i + 1 until variants.size) {
            val other = variants[j]
            if (other.minStart > variant.maxStart + maxSeekDistance) {
                break
            } else if (allFilters(other)) {
                result.addLast(other)
            }
        }

        return result
    }

}