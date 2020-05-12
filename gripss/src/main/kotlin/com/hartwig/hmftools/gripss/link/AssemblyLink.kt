package com.hartwig.hmftools.gripss.link

import com.hartwig.hmftools.gripss.StructuralVariantContext

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
                if (variants.size == 2 && variants[0].mateId?.equals(variants[1].vcfId) != true) {
                    for (variant in variants) {
                        result.compute(variant.vcfId) { _, v -> v?.plus(",$assembly") ?: assembly}
                    }
                }
            }
        }

        return result
    }

}


