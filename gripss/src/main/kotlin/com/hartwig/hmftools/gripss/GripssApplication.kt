package com.hartwig.hmftools.gripss

import htsjdk.variant.vcf.VCFFileReader
import java.io.File

fun main(args: Array<String>) {

    val time = System.currentTimeMillis();
    println("Starting")

    val inputVCF = "/Users/jon/hmf/analysis/gridss/CPCT02010893R_CPCT02010893T.gridss.vcf.gz"
//    val inputVCF = "/Users/jon/hmf/analysis/gridss/CPCT02010893T.gridss.somatic.full.vcf.gz"
    val outputVCF = "/Users/jon/hmf/analysis/gridss/CPCT02010893T.post.vcf"
    val filterConfig = GripssFilterConfig(0.03, 8, 0.005, 0.95, 350, 1000, 50)
    val config = GripssConfig(inputVCF, outputVCF, filterConfig)

    GripssApplication(config).use { x -> x.run() }

    println("Finished in ${(System.currentTimeMillis() - time) / 1000} seconds")
}


class GripssApplication(private val config: GripssConfig) : AutoCloseable, Runnable {

    private val fileReader = VCFFileReader(File(config.inputVcf), false)
    private val fileWriter = GripssVCF(config.outputVcf)

    override fun run() {
        fileWriter.writeHeader(fileReader.fileHeader)
        for (variantContext in fileReader) {

            val structuralVariant = StructuralVariantContext(variantContext)
            if (!structuralVariant.isHardFilter(config.filterConfig)) {
                fileWriter.writeVariant(structuralVariant.context(config.filterConfig))
            }
        }
    }

    override fun close() {
        fileReader.close()
        fileWriter.close()
    }
}

