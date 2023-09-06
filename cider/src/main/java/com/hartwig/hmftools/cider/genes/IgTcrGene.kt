package com.hartwig.hmftools.cider

import com.hartwig.hmftools.cider.genes.GenomicLocation

// see https://www.imgt.org/IMGTScientificChart/SequenceDescription/IMGTfunctionality.html
enum class IgTcrFunctionality
{
    FUNCTIONAL,
    ORF, // open reading frame
    PSEUDOGENE; // pseudogene

    fun toCode() : String
    {
        return when (this)
        {
            FUNCTIONAL -> "F"
            ORF -> "ORF"
            PSEUDOGENE -> "P"
        }
    }

    companion object
    {
        fun fromCode(code: String) : IgTcrFunctionality
        {
            return when (code)
            {
                "F" -> FUNCTIONAL
                "ORF" -> ORF
                "P" -> PSEUDOGENE
                else -> throw IllegalArgumentException("invalid IgTcrFunctionality code: $code")
            }
        }
    }
}

enum class IgTcrRegion
{
    V_REGION, D_REGION, J_REGION, CONSTANT
}

//
data class IgTcrGene(
    val geneName: String,
    val allele: String, // 01
    val region: IgTcrRegion,
    val functionality: IgTcrFunctionality,
    val geneLocation: GenomicLocation?,
    val anchorSequence: String?, // only valid for V / J gene
    val anchorLocation: GenomicLocation?
)
{
    val geneAllele: String get() { return "$geneName*$allele" }
    val isFunctional: Boolean get() { return functionality == IgTcrFunctionality.FUNCTIONAL }
}
