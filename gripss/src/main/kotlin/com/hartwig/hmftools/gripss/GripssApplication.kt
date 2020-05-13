package com.hartwig.hmftools.gripss

import com.hartwig.hmftools.gripss.link.LocalLink
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
        val structuralVariants: MutableList<StructuralVariantContext> = mutableListOf()
        for (variantContext in fileReader) {
            val structuralVariant = StructuralVariantContext(variantContext)
            if (!structuralVariant.isHardFilter(config.filterConfig)) {
                structuralVariants.add(structuralVariant)
            }
        }


        println("LINKING")

        val links = LocalLink.create(structuralVariants)

        println("WRITING")
        structuralVariants.forEach { x -> fileWriter.writeVariant(x.context(config.filterConfig, links.link(x.vcfId), links.link(x.mateId))) }

    }

    override fun close() {
        fileReader.close()
        fileWriter.close()
    }
}

