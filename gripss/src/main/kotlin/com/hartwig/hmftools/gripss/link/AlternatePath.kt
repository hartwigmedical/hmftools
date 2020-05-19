package com.hartwig.hmftools.gripss.link

import com.hartwig.hmftools.gripss.StructuralVariantContext
import com.hartwig.hmftools.gripss.store.LinkStore
import com.hartwig.hmftools.gripss.store.VariantStore
import org.apache.logging.log4j.LogManager
import java.util.*
import kotlin.collections.HashSet

data class AlternatePath(val vcfId: String, val path: List<Link>) {

    companion object {
        private val logger = LogManager.getLogger(this::class.java)

        operator fun invoke(assemblyLinkStore: LinkStore, variantStore: VariantStore, variants: List<StructuralVariantContext>): Collection<AlternatePath> {

            val failed = HashSet<String>()
            val result = HashMap<String, AlternatePath>()
            val transitiveLink = TransitiveLink(assemblyLinkStore, variantStore)

            for (variant in variants) {
                if (variant.mateId != null && !result.keys.contains(variant.vcfId) && !failed.contains(variant.vcfId)) {
                    val links = transitiveLink.transitiveLink(variant)
                    if (links.isNotEmpty()) {
                        val alternatePath = AlternatePath(variant.vcfId, links)
                        val reverseAlternatePath = AlternatePath(variant.mateId, links.map { x -> x.reverse() }.reversed())
                        result[variant.vcfId] = alternatePath
                        result[variant.mateId] = reverseAlternatePath
                        logger.info("Found alternate mapping of $variant CIPOS:${variant.confidenceInterval} IMPRECISE:${variant.imprecise} -> ${alternatePath.path()}")
                    } else {
                        failed.add(variant.mateId)
                    }
                }
            }

            return result.values
        }
    }


    fun transitiveLinks(): List<Link> {
        return path.filter { x -> x.link.startsWith("trs") }
    }

    fun path(): String {
        val stringJoiner = StringJoiner("")

        for (i in path.indices) {
            if (i == 0) {
                stringJoiner.add(path[i].vcfId)
            }

            val link = path[i]
            if (link.link == "PAIR") {
                stringJoiner.add("-")
            } else {
                stringJoiner.add("<${link.link}>")
            }

            stringJoiner.add(link.otherVcfId)
        }

        return stringJoiner.toString()
    }

}