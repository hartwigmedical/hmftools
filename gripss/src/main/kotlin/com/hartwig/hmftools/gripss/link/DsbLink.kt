package com.hartwig.hmftools.gripss.link

import com.hartwig.hmftools.gripss.StructuralVariantContext
import com.hartwig.hmftools.gripss.store.LinkStore
import com.hartwig.hmftools.gripss.store.VariantStore
import java.util.*
import kotlin.collections.HashMap
import kotlin.math.max

class DsbLink(private val duplicates: Set<String>, private val variants: List<StructuralVariantContext>, private val assemblyLinkStore: LinkStore) {

    companion object {
        const val MAX_DSB_DISTANCE = 30
        const val MAX_FIND_DISTANCE = 1000

        operator fun invoke(variantStore: VariantStore, assemblyLinkStore: LinkStore, duplicates: Set<String>): LinkStore {

            val linkByVariant = HashMap<String, Link>()
            val variants = variantStore.selectAll()
            val dsbLink = DsbLink(duplicates, variants, assemblyLinkStore)
            for (index in variants.indices) {
                val variant = variants[index]
                if (!linkByVariant.containsKey(variant.vcfId) && !duplicates.contains(variant.vcfId)) {
                    val links = dsbLink.dsbLinks(linkByVariant.size / 2 + 1, index)
                    links.forEach { x -> linkByVariant[x.vcfId] = x }
                }
            }
            return LinkStore(linkByVariant.values.associate { x -> Pair(x.vcfId, listOf(x)) })
        }
    }

    fun dsbLinks(linkId: Int, variantIndex: Int): List<Link> {
        val variant = variants[variantIndex]
        val nearby = findNearbyIndices(variantIndex)
        if (nearby.size == 1) {
            val otherIndex = nearby[0]
            val linkIsSymmetric = findNearbyIndices(otherIndex).size == 1
            if (!linkIsSymmetric) {
                return Collections.emptyList()
            }

            val other = variants[otherIndex]
            val alreadyLinked = assemblyLinkStore.linkedVariants(variant.vcfId).any { x -> x.otherVcfId == other.vcfId }
            if (alreadyLinked) {
                return Collections.emptyList()
            }

            return listOf(Link("dsb$linkId", Pair(variant, other)), Link("dsb$linkId", Pair(other, variant)))
        }

        return Collections.emptyList()
    }

    private fun findNearbyIndices(i: Int): List<Int> {
        val variant = variants[i]
        val result = mutableListOf<Int>()

        val idFilter = { other: StructuralVariantContext -> variant.vcfId != other.vcfId && variant.mateId?.equals(other.vcfId) != true }
        val overlapFilter = { other: StructuralVariantContext -> other.minStart <= variant.maxStart + MAX_DSB_DISTANCE && other.maxStart >= variant.minStart - MAX_DSB_DISTANCE }
        val orientationFilter = {other: StructuralVariantContext -> other.orientation != variant.orientation}
        val duplicateFilter = {other: StructuralVariantContext -> !duplicates.contains(other.vcfId)}
        val allFilters = {other: StructuralVariantContext -> idFilter(other) && overlapFilter(other) && orientationFilter(other) && duplicateFilter(other)}

        // Look forwards
        for (j in i + 1 until variants.size) {
            val other = variants[j]
            if (other.minStart > variant.maxStart + MAX_FIND_DISTANCE) {
                break
            } else if (allFilters(other)) {
                result.add(j)
            }
        }

        // Look backwards
        for (j in max(0, i - 1) downTo 0) {
            val other = variants[j]
            if (other.maxStart < variant.minStart - MAX_FIND_DISTANCE) {
                break
            } else if (allFilters(other)) {
                result.add(j)
            }
        }

        return result
    }
}