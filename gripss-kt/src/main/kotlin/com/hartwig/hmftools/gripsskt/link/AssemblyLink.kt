package com.hartwig.hmftools.gripsskt.link

import com.hartwig.hmftools.gripsskt.StructuralVariantContext
import com.hartwig.hmftools.gripsskt.store.LinkStore
import java.util.*
import kotlin.Comparator
import kotlin.collections.HashMap

class AssemblyLink {

    companion object {
        operator fun invoke(allVariants: Collection<StructuralVariantContext>): LinkStore {
            return LinkStore(AssemblyLink().links(allVariants))
        }
    }

    fun links(allVariants: Collection<StructuralVariantContext>): List<Link> {
        val variantsByAssembly: HashMap<String, MutableList<StructuralVariantContext>> = HashMap()
        for (variant in allVariants) {
            if (!variant.isSingle) {
                for (assembly in variant.assemblies()) {
                    variantsByAssembly.computeIfAbsent(assembly) { mutableListOf() }.add(variant)
                }
            }
        }

        val result = mutableListOf<Link>()
        for ((assembly, assemblyVariants) in variantsByAssembly) {
            if (assemblyVariants.size > 1) {
                result.addAll(createLinks(assembly, assemblyVariants))
            }
        }

        return result
    }

    private fun createLinks(assembly: String, variants: MutableList<StructuralVariantContext>): List<Link> {

        if (variants.size < 2) {
            return Collections.emptyList()
        }

        if (variants.size == 2) {
            return createLinks(assembly, Pair(variants[0], variants[1]))
        }

        var extraIdentifier = 1
        val result = mutableListOf<Link>()
        variants.sortWith(Comparator { x, y -> x.start.compareTo(y.start) })
        for (i in 0..variants.size - 2) {
            val current = variants[i]
            val next = variants[i + 1]

            val links = createLinks("${assembly}-${extraIdentifier}", Pair(current, next))
            if (links.isNotEmpty()) {
                extraIdentifier++
                result.addAll(links)
            }
        }

        return result
    }

    private fun createLinks(assembly: String, variants: Pair<StructuralVariantContext, StructuralVariantContext>): List<Link> {
        if (variants.first.mateId?.equals(variants.second.vcfId) != true) {
            return listOf(Link(assembly, Pair(variants.first, variants.second)), Link(assembly, Pair(variants.second, variants.first)))
        }
        return Collections.emptyList()
    }
}


