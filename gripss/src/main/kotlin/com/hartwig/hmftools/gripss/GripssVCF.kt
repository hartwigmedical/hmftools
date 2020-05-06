package com.hartwig.hmftools.gripss

import htsjdk.variant.variantcontext.VariantContext
import htsjdk.variant.variantcontext.writer.Options
import htsjdk.variant.variantcontext.writer.VariantContextWriterBuilder
import htsjdk.variant.vcf.VCFFilterHeaderLine
import htsjdk.variant.vcf.VCFHeader

const val MAX_NORMAL_SUPPORT = "maxNormalSupport";
const val MIN_NORMAL_COVERAGE = "minNormalCoverage";

class GripssVCF(outputVCF: String) : AutoCloseable {
    val writer = VariantContextWriterBuilder().setOutputFile(outputVCF).unsetOption(Options.INDEX_ON_THE_FLY).build()


    fun writeHeader(header: VCFHeader, config: GripssFilterConfig) {
        header.addMetaDataLine(VCFFilterHeaderLine(MAX_NORMAL_SUPPORT, "Filter variants with more than ${config.maxNormalSupport * 100}% of the supporting reads originating from the normal"))
        header.addMetaDataLine(VCFFilterHeaderLine(MIN_NORMAL_COVERAGE, "Filter variants with breakend coverage of less than ${config.minNormalCoverage} fragments coverage"))
        writer.writeHeader(header)
    }

    fun writeVariant(context: VariantContext) {
        writer.add(context)
    }

    override fun close() {
        println("Closing file")
        writer.close()
    }

}