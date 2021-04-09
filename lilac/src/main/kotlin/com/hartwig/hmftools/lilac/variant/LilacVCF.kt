package com.hartwig.hmftools.lilac.variant

import com.hartwig.hmftools.lilac.LilacApplication
import htsjdk.variant.variantcontext.VariantContext
import htsjdk.variant.variantcontext.writer.VariantContextWriterBuilder
import htsjdk.variant.vcf.*
import java.io.File

const val HLA = "HLA"


class LilacVCF(outputVCF: String, templateVCF: String) : AutoCloseable {
    companion object {
        const val HLA = "HLA"
    }

    private val header = VCFFileReader(File(templateVCF), false).use { it.fileHeader }
    private val writer = VariantContextWriterBuilder()
            .setReferenceDictionary(header.sequenceDictionary)
            .setOutputFile(outputVCF)
            .build()

    fun writeHeader() {
        val newHeader = VCFHeader(header)
        newHeader.addMetaDataLine(VCFHeaderLine("lilacVersion", LilacApplication.VERSION.version()))
        newHeader.addMetaDataLine(VCFInfoHeaderLine(HLA, VCFHeaderLineCount.UNBOUNDED, VCFHeaderLineType.String, "HLA Type"))

        writer.writeHeader(newHeader)
    }

    fun writeVariant(context: VariantContext) {
        writer.add(context)
    }

    override fun close() {
        writer.close()
    }

}