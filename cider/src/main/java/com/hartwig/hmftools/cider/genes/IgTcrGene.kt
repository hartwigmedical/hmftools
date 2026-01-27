package com.hartwig.hmftools.cider

import com.hartwig.hmftools.cider.genes.GenomicLocation
import com.hartwig.hmftools.cider.genes.GenomicLocation.Companion.toChrBaseRegion

// TODO: should get rid of this

// see https://www.imgt.org/IMGTScientificChart/SequenceDescription/IMGTfunctionality.html
typealias IgTcrFunctionality = com.hartwig.hmftools.common.cider.IgTcrFunctionality

typealias IgTcrRegion = com.hartwig.hmftools.common.cider.IgTcrRegion

// use kotlin data class to use the nicer kotlin functionality
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
    
    companion object {

        fun fromCommonIgTcrGene(commonGene: com.hartwig.hmftools.common.cider.IgTcrGene): IgTcrGene = IgTcrGene(
            geneName = commonGene.geneName,
            allele = commonGene.allele,
            region = commonGene.region,
            functionality = commonGene.functionality,
            geneLocation = GenomicLocation.fromChrBaseRegionStrand(commonGene.geneLocation, commonGene.geneStrand),
            anchorSequence = commonGene.anchorSequence,
            anchorLocation = GenomicLocation.fromChrBaseRegionStrand(commonGene.anchorLocation, commonGene.geneStrand)
        )

        fun toCommonIgTcrGene(gene: IgTcrGene): com.hartwig.hmftools.common.cider.IgTcrGene
        {
            return com.hartwig.hmftools.common.cider.IgTcrGene(
                gene.geneName,
                gene.allele,
                gene.region,
                gene.functionality,
                toChrBaseRegion(gene.geneLocation),
                gene.geneLocation?.strand,
                gene.geneLocation?.altAssemblyName,
                gene.anchorSequence,
                toChrBaseRegion(gene.anchorLocation)
            )
        }
    }
}
