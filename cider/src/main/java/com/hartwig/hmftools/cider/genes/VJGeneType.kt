package com.hartwig.hmftools.cider.genes

// this file defines the gene segments that encode the variable region
// of Ig/TCR locus. They are V, D, J and constant

// the various Ig/TCR locus
enum class IgTcrLocus
{
    IGH,
    IGK,
    IGL,
    TRA_TRD,
    TRB,
    TRG;

    fun prettyPrint() : String
    {
        return if (this == TRA_TRD)
            "TRA/TRD"
        else
            toString()
    }

    companion object
    {
        fun fromGeneName(geneName: String) : IgTcrLocus
        {
            val locusName = geneName.substring(0, 3)

            if (locusName == "TRA" || locusName == "TRD")
                return TRA_TRD

            return valueOf(locusName)
        }
    }
}

enum class VJ
{
    V, J
}

// make this a inner type of VJGene
// create a new type called IgLocus to house IGH, TRA, TRB etc
// note: IGH, TRA etc are locus, IGHV, TRAJ etc are gene segments
// this also include KDE, for simplicity we treat it just like a J anchor
//
enum class VJGeneType(val locus: IgTcrLocus, val vj: VJ)
{
    IGHV(IgTcrLocus.IGH, VJ.V),
    IGHJ(IgTcrLocus.IGH, VJ.J),
    IGKV(IgTcrLocus.IGK, VJ.V),
    IGKJ(IgTcrLocus.IGK, VJ.J),
    IGLV(IgTcrLocus.IGL, VJ.V),
    IGLJ(IgTcrLocus.IGL, VJ.J),
    TRAV(IgTcrLocus.TRA_TRD, VJ.V),
    TRAJ(IgTcrLocus.TRA_TRD, VJ.J),
    TRBV(IgTcrLocus.TRB, VJ.V),
    TRBJ(IgTcrLocus.TRB, VJ.J),
    TRDV(IgTcrLocus.TRA_TRD, VJ.V),
    TRDJ(IgTcrLocus.TRA_TRD, VJ.J),
    TRGV(IgTcrLocus.TRG, VJ.V),
    TRGJ(IgTcrLocus.TRG, VJ.J);

    // gene types that can be paired with
    fun pairedVjGeneTypes() : List<VJGeneType>
    {
        return when (this)
        {
            IGHV -> listOf(IGHJ)
            IGHJ -> listOf(IGHV)
            IGKV -> listOf(IGKJ)
            IGKJ -> listOf(IGKV)
            IGLV -> listOf(IGLJ)
            IGLJ -> listOf(IGLV)
            TRAV -> listOf(TRAJ, TRDJ)
            TRAJ -> listOf(TRAV, TRDV)
            TRBV -> listOf(TRBJ)
            TRBJ -> listOf(TRBV)
            TRDV -> listOf(TRDJ, TRAJ)
            TRDJ -> listOf(TRDV, TRAV)
            TRGV -> listOf(TRGJ)
            TRGJ -> listOf(TRGV)
        }
    }

    companion object
    {
        val IGKINTR = "IGKINTR"
        val IGKDEL = "IGKDEL"

        fun fromGeneName(geneName: String) : VJGeneType
        {
            return when (geneName)
            {
                IGKINTR -> IGKV
                IGKDEL -> IGKJ
                else -> valueOf(geneName.substring(0, 4))
            }
        }
    }
}
