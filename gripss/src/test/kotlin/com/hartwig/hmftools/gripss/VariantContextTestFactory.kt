package com.hartwig.hmftools.gripss

import com.hartwig.hmftools.extensions.CIPOS
import com.hartwig.hmftools.extensions.CIRPOS
import htsjdk.variant.variantcontext.GenotypeBuilder
import htsjdk.variant.variantcontext.VariantContext
import htsjdk.variant.variantcontext.VariantContextBuilder
import htsjdk.variant.vcf.VCFCodec
import htsjdk.variant.vcf.VCFHeader
import htsjdk.variant.vcf.VCFHeaderVersion
import kotlin.math.abs

object VariantContextTestFactory {

    private val codec: VCFCodec by lazy {
        val codec = VCFCodec()
        val header = VCFHeader(setOf(), setOf("NORMAL", "TUMOR"))
        codec.setVCFHeader(header, VCFHeaderVersion.VCF4_2)
        return@lazy codec
    }

    fun decode(line: String): VariantContext {
        return codec.decode(line)
    }

    fun createVariant(contig: String, pos: Int, id: String, ref: String, alt: String, qual: Int, filter: Collection<String>): VariantContext {
        val line = "$contig\t${pos}\t$id\t${ref}\t${alt}\t${qual}\t${filter.joinToString(",")}\tinfo\tGT:BVF:VF:REF:REFPAIR\t./.:1:1:1:1\t./.:10:10:1:1"
        return decode(line)
    }

    fun createVariant(contig: String, pos: Int, id: String, ref: String, alt: String, qual: Int, mateId: String): VariantContext {
        val line = "$contig\t${pos}\t$id\t${ref}\t${alt}\t${qual}\t.\tMATEID=$mateId\tGT:BVF:VF:REF:REFPAIR\t./.:1:1:1:1\t./.:10:10:1:1"
        return decode(line)
    }

    fun createImpreciseVariant(contig: String, pos: Int, id: String, ref: String, alt: String, qual: Int, mateId: String): VariantContext {
        val line = "$contig\t${pos}\t$id\t${ref}\t${alt}\t${qual}\t.\tIMPRECISE;MATEID=$mateId\tGT:BVF:VF:REF:REFPAIR\t./.:1:1:1:1\t./.:10:10:1:1"
        return decode(line)
    }

    fun VariantContext.toSv(): StructuralVariantContext = StructuralVariantContext(this)

    fun VariantContext.cipos(cipos: Pair<Int, Int>, cirpos: Pair<Int, Int>): VariantContext {
        val builder = VariantContextBuilder(this)
        builder.attribute(CIPOS, listOf(-abs(cipos.first), abs(cipos.second)))
        builder.attribute(CIRPOS, listOf(-abs(cirpos.first), abs(cirpos.second)))
        return builder.make()
    }

    fun VariantContext.setAttribute(attribute: String, value: Any): VariantContext {
        val builder = VariantContextBuilder(this)
        builder.attribute(attribute, value)
        return builder.make()
    }

    fun VariantContext.addGenotypeAttribute(attribute: String, normal: Int, tumor: Int): VariantContext {
        val normalBuilder = GenotypeBuilder(this.getGenotype(0))
        normalBuilder.attribute(attribute, normal)

        val tumorBuilder = GenotypeBuilder(this.getGenotype(1))
        tumorBuilder.attribute(attribute, tumor)

        val builder = VariantContextBuilder(this)
        builder.genotypes(listOf(normalBuilder.make(), tumorBuilder.make()))
        return builder.make()
    }

    fun VariantContext.splitReads(normal: Int, tumor: Int): VariantContext {
        return addGenotypeAttribute("SR", normal, tumor)
    }

    fun VariantContext.fragmentSupport(normal: Int, tumor: Int): VariantContext {
        return addGenotypeAttribute("BVF", normal, tumor).addGenotypeAttribute("VF", normal, tumor)
    }
}

