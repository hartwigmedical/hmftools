package com.hartwig.hmftools.gripss

import com.hartwig.hmftools.extensions.dedup.AssemblyDedup
import com.hartwig.hmftools.gripss.link.AssemblyLink
import com.hartwig.hmftools.gripss.store.LinkStore
import com.hartwig.hmftools.gripss.store.VariantStore
import htsjdk.variant.vcf.VCFFileReader
import java.io.File

fun main(args: Array<String>) {

    val time = System.currentTimeMillis();
    println("Starting")

    val inputVCF = "/Users/jon/hmf/analysis/gridss/CPCT02010893R_CPCT02010893T.gridss.vcf.gz"
//    val inputVCF = "/Users/jon/hmf/analysis/gridss/CPCT02010893T.gridss.somatic.vcf"
    val outputVCF = "/Users/jon/hmf/analysis/gridss/CPCT02010893T.post.vcf"
    val filterConfig = GripssFilterConfig(
            0.03,
            8,
            0.005,
            0.95,
            1000,
            350,
            50,
            6,
            50,
            5,
            32)
    val config = GripssConfig(inputVCF, outputVCF, filterConfig)

    GripssApplication(config).use { x -> x.run() }

    println("Finished in ${(System.currentTimeMillis() - time) / 1000} seconds")
}


class GripssApplication(private val config: GripssConfig) : AutoCloseable, Runnable {

    private val fileReader = VCFFileReader(File(config.inputVcf), false)
    private val fileWriter = GripssVCF(config.outputVcf)

    override fun run() {
        println("READING")

        fileWriter.writeHeader(fileReader.fileHeader)
        val structuralVariants = hardFilterVariants(fileReader)


        println("LINKING")

        val variantStore = VariantStore.create(structuralVariants)
        val assemblyLinks = AssemblyLink().create(structuralVariants)
        val links = LinkStore.create(assemblyLinks)
        val assemblyDedup = AssemblyDedup(links, variantStore);

        for (variant in structuralVariants) {
            assemblyDedup.dedup(variant)
        }

        println("WRITING")
        structuralVariants.forEach { x -> fileWriter.writeVariant(x.context(config.filterConfig, links.localLinkedBy(x.vcfId), links.localLinkedBy(x.mateId))) }

    }

    fun hardFilterVariants(fileReader: VCFFileReader): List<StructuralVariantContext> {
        val unfiltered: MutableSet<String> = mutableSetOf()
        val hardFilter: MutableSet<String> = mutableSetOf()
        val structuralVariants: MutableList<StructuralVariantContext> = mutableListOf()

        for (variantContext in fileReader) {
            val structuralVariant = StructuralVariantContext(variantContext)
            if (hardFilter.contains(structuralVariant.vcfId) || structuralVariant.isHardFilter(config.filterConfig)) {
                structuralVariant.mateId?.let { hardFilter.add(it) }
            } else {
                unfiltered.add(variantContext.id)
                structuralVariants.add(structuralVariant)
            }
        }

        val mateIsValidOrNull = { x: StructuralVariantContext -> x.mateId?.let { unfiltered.contains(it) } != false}
        return structuralVariants.filter { x -> !hardFilter.contains(x.vcfId) && mateIsValidOrNull(x) }
    }


    override fun close() {
        fileReader.close()
        fileWriter.close()
    }
}

