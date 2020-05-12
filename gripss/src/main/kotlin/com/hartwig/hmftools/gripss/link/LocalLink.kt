package com.hartwig.hmftools.gripss.link

import com.hartwig.hmftools.gripss.StructuralVariantContext
import org.apache.logging.log4j.util.Strings

data class LocalLink(val linkById: Map<String, String>) {

    companion object Factory {

        fun create(variants: Collection<StructuralVariantContext>): LocalLink {

            val assemblyLink = AssemblyLink()
            for (variant in variants) {
                assemblyLink.addVariant(variant)
            }

            val assemblyLinkMap = assemblyLink.createMap()

            return LocalLink(assemblyLinkMap)
        }

    }

    fun link(id: String?): String {
        return linkById.getOrDefault(id, Strings.EMPTY)
    }


}