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

    fun selectOthersNearby(variant: StructuralVariantContext, maxDistance: Pair<Int, Int>, filter: (StructuralVariantContext) -> Boolean = { _ -> true }): List<StructuralVariantContext> {
        val minStart = variant.minStart - maxDistance.first
        val maxStart = variant.maxStart + maxDistance.second
        val idFilter = { other: StructuralVariantContext -> variant.vcfId != other.vcfId && variant.mateId?.equals(other.vcfId) != true }
        val overlapFilter = { other: StructuralVariantContext -> other.minStart <= maxStart && other.maxStart >= minStart }
        return variantsByChromosome.getOrDefault(variant.contig, Collections.emptyList()).filter { x -> idFilter(x) && overlapFilter(x) && filter(x) }
    }

    fun selectOthersNearby(variant: StructuralVariantContext, maxDistance: Int, filter: (StructuralVariantContext) -> Boolean = { _ -> true }): List<StructuralVariantContext> {
        return selectOthersNearby(variant, Pair(maxDistance, maxDistance), filter)
    }

}