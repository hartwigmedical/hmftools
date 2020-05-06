package com.hartwig.hmftools.gripss

import htsjdk.variant.variantcontext.VariantContext
import htsjdk.variant.variantcontext.writer.Options
import htsjdk.variant.variantcontext.writer.VariantContextWriterBuilder
import htsjdk.variant.vcf.VCFFilterHeaderLine
import htsjdk.variant.vcf.VCFHeader
import htsjdk.variant.vcf.VCFHeaderLineType
import htsjdk.variant.vcf.VCFInfoHeaderLine

const val MAX_NORMAL_SUPPORT = "maxNormalSupport";
const val MIN_NORMAL_COVERAGE = "minNormalCoverage";
const val MIN_TUMOR_AF = "minTumorAF";
const val SHORT_STRAND_BIAS = "shortStrandBias";
const val TAF = "TAF";

class GripssVCF(outputVCF: String) : AutoCloseable {
    val writer = VariantContextWriterBuilder().setOutputFile(outputVCF).unsetOption(Options.INDEX_ON_THE_FLY).build()


    fun writeHeader(header: VCFHeader, config: GripssFilterConfig) {
        header.addMetaDataLine(VCFFilterHeaderLine(MIN_TUMOR_AF, "Filter variants with less than ${config.minTumorAF * 100}% tumor allelic frequency"))
        header.addMetaDataLine(VCFFilterHeaderLine(MAX_NORMAL_SUPPORT, "Filter variants with more than ${config.maxNormalSupport * 100}% of the supporting reads originating from the normal"))
        header.addMetaDataLine(VCFFilterHeaderLine(MIN_NORMAL_COVERAGE, "Filter variants with breakend coverage of less than ${config.minNormalCoverage} fragments coverage"))
        header.addMetaDataLine(VCFFilterHeaderLine(SHORT_STRAND_BIAS, "Filter deletion or duplication breakpoints under 1000bp with a split read strand bias of more than ${config.maxShortStrandBias}"))
        header.addMetaDataLine(VCFInfoHeaderLine(TAF, 1, VCFHeaderLineType.Float, "Description"))
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