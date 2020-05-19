package com.hartwig.hmftools.gripss.store

import com.hartwig.hmftools.gripss.link.Link
import org.apache.logging.log4j.util.Strings
import java.util.*
import kotlin.collections.HashMap

class LinkStore(private val variantsByLink: Map<String, List<String>>, private val linkByVariant: Map<String, List<String>>) {

    companion object Factory {
        fun create(links: Collection<Link>): LinkStore {
            val variantsByLink = HashMap<String, MutableList<String>>()
            val linkByVariant = HashMap<String, MutableList<String>>()

            for ((link, vcfId) in links) {
                linkByVariant.computeIfAbsent(vcfId) { mutableListOf()}.add(link)
                variantsByLink.computeIfAbsent(link) { mutableListOf()}.add(vcfId)
            }
            return LinkStore(variantsByLink, linkByVariant)
        }
    }

    fun localLinkedBy(vcfId: String?): String {
        return linkByVariant.get(vcfId)?.joinToString(",")  ?: Strings.EMPTY
    }

    fun linkedVariants(vcfId: String): List<String> {
        return linksByVariant(vcfId).flatMap { x -> variantsByLink(x) }.filter { x -> x != vcfId }
    }

    private fun variantsByLink(link: String): List<String> {
        return variantsByLink.getOrDefault(link, Collections.emptyList())
    }

    private fun linksByVariant(vcfId: String): List<String> {
        return linkByVariant.getOrDefault(vcfId, Collections.emptyList())
    }

}