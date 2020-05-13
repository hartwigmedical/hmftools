package com.hartwig.hmftools.gripss.link

import com.hartwig.hmftools.gripss.StructuralVariantContext
import java.util.*
import kotlin.Comparator
import kotlin.collections.HashMap

class AssemblyLink {

    private val assemblyMap: HashMap<String, MutableList<StructuralVariantContext>> = HashMap()

    fun addVariant(variant: StructuralVariantContext) {
        if (!variant.isSingle) {
            for (assembly in variant.assemblies()) {
                assemblyMap.computeIfAbsent(assembly) { mutableListOf() }.add(variant)
            }
        }
    }

    fun createMap(): Map<String, String> {
        val result: HashMap<String, String> = HashMap()

        for ((assembly, variants) in assemblyMap) {
            if (variants.size > 1) {
                result.mergeAssemblies(createMap(assembly, variants))
            }
        }

        return result
    }

    private fun createMap(assembly: String, variants: MutableList<StructuralVariantContext>): Map<String, String> {
        if (variants.size < 2) {
            return Collections.emptyMap()
        }

        if (variants.size == 2) {
            return createMap(assembly, Pair(variants[0], variants[1]));
        }

        var extraIdentifier = 1
        val result: HashMap<String, String> = HashMap()
        variants.sortWith(Comparator { x, y -> x.start.compareTo(y.start) })
        for (i in 0..variants.size - 2) {
            val current = variants[i]
            val next = variants[i + 1]

            val pairMap = createMap("${assembly}-${extraIdentifier}", Pair(current, next))
            if (pairMap.isNotEmpty()) {
                extraIdentifier++;
                result.mergeAssemblies(pairMap)
            }
        }

        return result
    }

    private fun createMap(assembly: String, variants: Pair<StructuralVariantContext, StructuralVariantContext>): Map<String, String> {
        val result: HashMap<String, String> = HashMap()
        if (variants.first.mateId?.equals(variants.second.vcfId) != true) {
            result[variants.first.vcfId] = assembly
            result[variants.second.vcfId] = assembly
        }
        return result
    }

    private fun HashMap<String,String>.mergeAssemblies(other: Map<String, String>) {
        for ((vcfId, assembly) in other) {
            this.compute(vcfId) { _, v -> v?.plus(",$assembly") ?: assembly }
        }
    }

}


