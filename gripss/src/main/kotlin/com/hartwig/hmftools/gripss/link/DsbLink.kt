package com.hartwig.hmftools.gripss.link

import com.hartwig.hmftools.gripss.StructuralVariantContext
import com.hartwig.hmftools.gripss.store.LinkStore
import com.hartwig.hmftools.gripss.store.VariantStore
import java.util.*
import kotlin.collections.HashMap

class DsbLink(private val variantStore: VariantStore, private val assemblyLinkStore: LinkStore, private val duplicates: Set<String>) {

    companion object {
        const val MAX_DSB_DISTANCE = 30

        operator fun invoke(variantStore: VariantStore, assemblyLinkStore: LinkStore, duplicates: Set<String>, variants: List<StructuralVariantContext>): LinkStore {

            val linkByVariant = HashMap<String, Link>()
            val dsbLink = DsbLink(variantStore, assemblyLinkStore, duplicates)
            for (variant in variants) {
                if (!linkByVariant.containsKey(variant.vcfId) && !duplicates.contains(variant.vcfId)) {
                    val links = dsbLink.dsbLinks(linkByVariant.size / 2 + 1, variant)
                    links.forEach { x -> linkByVariant[x.vcfId] = x }
                }
            }
            return LinkStore(linkByVariant.values.associate { x -> Pair(x.vcfId, listOf(x)) })
        }
    }

    fun dsbLinks(i: Int, variant: StructuralVariantContext): List<Link> {
        val nearby = variantStore.selectDoubleStrandedBreaks(variant)
        if (nearby.size == 1) {
            val other = nearby[0]
            val linkIsSymmetric = variantStore.selectDoubleStrandedBreaks(other).size == 1
            if (!linkIsSymmetric) {
                return Collections.emptyList()
            }

            val alreadyLinked = assemblyLinkStore.linkedVariants(variant.vcfId).filter { x -> x.otherVcfId == other.vcfId }.isNotEmpty()
            if (alreadyLinked) {
                return Collections.emptyList()
            }

            return listOf(Link("dsb$i", Pair(variant, other)), Link("dsb$i", Pair(other, variant)))
        }

        return Collections.emptyList()
    }

    private fun VariantStore.selectDoubleStrandedBreaks(variant: StructuralVariantContext): List<StructuralVariantContext> {
        val dsbFilter = { other: StructuralVariantContext -> other.orientation != variant.orientation && !duplicates.contains(other.vcfId) }
        return selectOthersNearby(variant, MAX_DSB_DISTANCE, dsbFilter)
    }

}