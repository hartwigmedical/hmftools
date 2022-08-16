package com.hartwig.hmftools.cdr3

data class IgConstantRegion(val type: Type, val geneLocation: GeneLocation)
{
    enum class Type
    {
        IGHA, IGHG, IGHD, IGHE, IGHM,
        IGKC, IGLC, TRAC, TRBC, TRDC, TRGC
    }

    fun getCorrespondingJ() : VJGeneType
    {
        return when (type)
        {
            Type.IGHA, Type.IGHG, Type.IGHD, Type.IGHE, Type.IGHM -> VJGeneType.IGHJ
            Type.IGKC -> VJGeneType.IGKJ
            Type.IGLC -> VJGeneType.IGLJ
            Type.TRAC -> VJGeneType.TRAJ
            Type.TRBC -> VJGeneType.TRBJ
            Type.TRDC -> VJGeneType.TRDJ
            Type.TRGC -> VJGeneType.TRGJ
        }
    }
}
