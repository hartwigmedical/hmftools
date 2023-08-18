package com.hartwig.hmftools.cider.genes

import com.hartwig.hmftools.cider.IgTcrLocus
import com.hartwig.hmftools.cider.VJGeneType

data class IgTcrConstantRegion(val locus: IgTcrLocus, val genomeLocation: GenomicLocation)
{
    fun getCorrespondingJ() : List<VJGeneType>
    {
        return when (locus)
        {
            IgTcrLocus.IGH -> listOf(VJGeneType.IGHJ)
            IgTcrLocus.IGK -> listOf(VJGeneType.IGKJ)
            IgTcrLocus.IGL -> listOf(VJGeneType.IGLJ)
            IgTcrLocus.TRB -> listOf(VJGeneType.TRBJ)
            IgTcrLocus.TRG -> listOf(VJGeneType.TRGJ)
            IgTcrLocus.TRA_TRD -> listOf(VJGeneType.TRAJ, VJGeneType.TRDJ)
        }
    }
}
