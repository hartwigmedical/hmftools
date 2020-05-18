package com.hartwig.hmftools.gripss.store

import com.hartwig.hmftools.gripss.StructuralVariantContext
import java.util.*
import kotlin.collections.HashMap

class VariantStore(
        private val variants: List<StructuralVariantContext>,
        private val variantsById: Map<String, StructuralVariantContext>,
        private val variantsByChromosome: Map<String, List<StructuralVariantContext>>) {

    companion object Factory {
        fun create(variants: List<StructuralVariantContext>): VariantStore {
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

    fun selectNearbyFacingVariants(variant: StructuralVariantContext, maxDistance: Int = 1000): List<StructuralVariantContext> {
        val commonFilter: (StructuralVariantContext) -> Boolean = { x -> x.orientation != variant.orientation}
        val filter: (StructuralVariantContext) -> Boolean = if (variant.orientation == 1.toByte())
            { x -> commonFilter(x) && x.start <= variant.start && x.start >= variant.start - maxDistance}
        else
            { x -> commonFilter(x) && x.start >= variant.start && x.start <= variant.start + maxDistance}
        return variantsByChromosome.getOrDefault(variant.contig, Collections.emptyList()).filter(filter)
    }

    fun selectAlternatives(variant: StructuralVariantContext): List<StructuralVariantContext> {
        return selectInRange(variant) { x -> x.vcfId != variant.vcfId && x.mateId?.equals(variant.vcfId) != true && x.orientation == variant.orientation }
    }

    private fun selectInRange(variant: StructuralVariantContext, filter: (StructuralVariantContext) -> Boolean): List<StructuralVariantContext> {
        return variantsByChromosome.getOrDefault(variant.contig, Collections.emptyList()).filter { other ->
            return@filter filter(other) && variant.confidenceIntervalsOverlap(other)
        }
    }


}