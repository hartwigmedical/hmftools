package com.hartwig.hmftools.lilac.variant

import com.hartwig.hmftools.lilac.LilacApplication
import com.hartwig.hmftools.lilac.hla.HlaAllele
import htsjdk.variant.variantcontext.VariantContext
import htsjdk.variant.variantcontext.VariantContextBuilder
import htsjdk.variant.variantcontext.writer.VariantContextWriterBuilder
import htsjdk.variant.vcf.*
import java.io.File

class LilacVCF(outputVCF: String, templateVCF: String) : AutoCloseable {
    companion object {
        const val HLA = "HLA"
    }

    private val header = VCFFileReader(File(templateVCF), false).use { it.fileHeader }
    private val writer = VariantContextWriterBuilder()
            .setReferenceDictionary(header.sequenceDictionary)
            .setOutputFile(outputVCF)
            .build()

    fun writeHeader(): LilacVCF {
        val newHeader = VCFHeader(header)
        newHeader.addMetaDataLine(VCFHeaderLine("lilacVersion", LilacApplication.VERSION.version()))
        newHeader.addMetaDataLine(VCFInfoHeaderLine(HLA, VCFHeaderLineCount.UNBOUNDED, VCFHeaderLineType.String, "HLA Type"))

        writer.writeHeader(newHeader)
        return this
    }

    fun writeVariant(context: VariantContext, alleles: Set<HlaAllele>) {
        val hlaLabel = if (alleles.isEmpty()) {
            listOf("UNKNOWN")
        } else {
            alleles.map { it.toString() }
        }

        val newContext = VariantContextBuilder(context).attribute(HLA, hlaLabel).make()
        writer.add(newContext)
    }

    override fun close() {
        writer.close()
    }

}