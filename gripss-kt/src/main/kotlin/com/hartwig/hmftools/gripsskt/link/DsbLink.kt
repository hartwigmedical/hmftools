package com.hartwig.hmftools.gripsskt.link

import com.hartwig.hmftools.gripsskt.StructuralVariantContext
import com.hartwig.hmftools.gripsskt.store.LinkStore
import com.hartwig.hmftools.gripsskt.store.VariantStore
import java.util.*
import kotlin.collections.HashMap

class DsbLink(private val duplicates: Set<String>, private val variantStore: VariantStore, private val assemblyLinkStore: LinkStore) {

    companion object {
        private const val MAX_DSB_SEEK_DISTANCE = 1000
        private const val MAX_DSB_ADDITIONAL_DISTANCE = 30

        operator fun invoke(variantStore: VariantStore, assemblyLinkStore: LinkStore, duplicates: Set<String>): LinkStore {

            val linkByVariant = HashMap<String, Link>()
            val variants = variantStore.selectAll()
            val dsbLink = DsbLink(duplicates, variantStore, assemblyLinkStore)
            for (variant in variants) {
                if (!linkByVariant.containsKey(variant.vcfId) && !duplicates.contains(variant.vcfId)) {
                    val links = dsbLink.dsbLinks(linkByVariant.size / 2 + 1, variant)
                    links.forEach { x -> linkByVariant[x.vcfId] = x }
                }
            }

            return LinkStore(linkByVariant.values.associate { x -> Pair(x.vcfId, listOf(x)) })
        }
    }

    fun dsbLinks(linkId: Int, variant: StructuralVariantContext): List<Link> {
        val nearby = findNearby(variant)
        if (nearby.size == 1) {
            val other = nearby.iterator().next()

            val linkIsSymmetric = findNearby(other).size == 1
            if (!linkIsSymmetric) {
                return Collections.emptyList()
            }

            val alreadyLinked = assemblyLinkStore.linkedVariants(variant.vcfId).any { x -> x.otherVcfId == other.vcfId }
            if (alreadyLinked) {
                return Collections.emptyList()
            }

            return listOf(Link("dsb$linkId", Pair(variant, other)), Link("dsb$linkId", Pair(other, variant)))
        }

        return Collections.emptyList()
    }

    private fun findNearby(variant: StructuralVariantContext): Collection<StructuralVariantContext> {
        val orientationFilter = { other: StructuralVariantContext -> other.orientation != variant.orientation }
        val duplicateFilter = { other: StructuralVariantContext -> !duplicates.contains(other.vcfId) }
        val allFilters = { other: StructuralVariantContext -> orientationFilter(other) && duplicateFilter(other) }
        return variantStore.selectOthersNearby(variant, MAX_DSB_ADDITIONAL_DISTANCE, MAX_DSB_SEEK_DISTANCE, allFilters)
    }
}