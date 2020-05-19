package com.hartwig.hmftools.gripss.store

import com.hartwig.hmftools.gripss.link.Link
import org.apache.logging.log4j.util.Strings
import java.util.*
import kotlin.collections.HashMap

class LinkStore(private val linkByVariant: Map<String, List<Link>>) {

    companion object Factory {
        fun create(links: Collection<Link>): LinkStore {
            val variantsByLink = HashMap<String, MutableList<Link>>()
            val linkByVariant = HashMap<String, MutableList<Link>>()

            for (link in links) {
                linkByVariant.computeIfAbsent(link.vcfId) { mutableListOf()}.add(link)
                variantsByLink.computeIfAbsent(link.link) { mutableListOf()}.add(link)
            }
            return LinkStore(linkByVariant)
        }
    }

    fun localLinkedBy(vcfId: String?): String {
        return linkByVariant[vcfId]?.joinToString(",") { x -> x.link } ?: Strings.EMPTY
    }

    fun linkedVariants(vcfId: String): List<Link> {
        return linkByVariant.getOrDefault(vcfId, Collections.emptyList())
    }
}