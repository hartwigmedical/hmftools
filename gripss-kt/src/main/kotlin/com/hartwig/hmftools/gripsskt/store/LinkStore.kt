package com.hartwig.hmftools.gripsskt.store

import com.hartwig.hmftools.gripsskt.link.Link
import org.apache.logging.log4j.util.Strings
import java.util.*
import kotlin.collections.HashMap

class LinkStore(private val linkByVariant: Map<String, List<Link>>) {

    companion object Factory {
        operator fun invoke(links: Collection<Link>): LinkStore {
            val linkByVariant = HashMap<String, MutableList<Link>>()

            for (link in links) {
                linkByVariant.computeIfAbsent(link.vcfId) { mutableListOf() }.add(link)
            }
            return LinkStore(linkByVariant)
        }

        operator fun invoke(vararg stores: LinkStore): LinkStore {
            val linkByVariant = HashMap<String, MutableList<Link>>()

            for (store in stores) {
                for (link in store.linkByVariant.values.flatten()) {
                    linkByVariant.computeIfAbsent(link.vcfId) { mutableListOf() }.add(link)
                }
            }

            return LinkStore(linkByVariant)
        }
    }

    operator fun get(vcfId: String?): String {
        return linkByVariant[vcfId]?.joinToString(",") { x -> x.link } ?: Strings.EMPTY
    }

    fun linkedVariants(vcfId: String): List<Link> {
        return linkByVariant.getOrDefault(vcfId, Collections.emptyList())
    }

    fun linkedVariants(): Set<String> {
        return linkByVariant.keys
    }
}