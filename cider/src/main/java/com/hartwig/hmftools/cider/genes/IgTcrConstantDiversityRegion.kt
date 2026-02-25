package com.hartwig.hmftools.cider.genes

data class IgTcrConstantDiversityRegion(val locus: IgTcrLocus, val genomeLocation: GenomicLocation, val geneName: String)
{
    init
    {
        require(genomeLocation.inPrimaryAssembly)
    }

    fun getCorrespondingVJ() : List<VJGeneType>
    {
        return when (locus)
        {
            IgTcrLocus.IGH -> listOf(VJGeneType.IGHV, VJGeneType.IGHJ)
            IgTcrLocus.IGK -> listOf(VJGeneType.IGKV, VJGeneType.IGKJ)
            IgTcrLocus.IGL -> listOf(VJGeneType.IGLV, VJGeneType.IGLJ)
            IgTcrLocus.TRB -> listOf(VJGeneType.TRBV, VJGeneType.TRBJ)
            IgTcrLocus.TRG -> listOf(VJGeneType.TRGV, VJGeneType.TRGJ)
            IgTcrLocus.TRA_TRD -> listOf(VJGeneType.TRAV, VJGeneType.TRDV, VJGeneType.TRAJ, VJGeneType.TRDJ)
        }
    }
}
