package com.hartwig.hmftools.cider

// this file defines the gene segments that encode the variable region
// of Ig/TCR locus. They are V, D, J and constant

// the various Ig/TCR locus
enum class IgTcrLocus
{
    IGH,
    IGK,
    IGL,
    TRA,
    TRB,
    TRD,
    TRG
}

enum class VJ
{
    V, J
}

// make this a inner type of VJGene
// create a new type called IgLocus to house IGH, TRA, TRB etc
// note: IGH, TRA etc are locus, IGHV, TRAJ etc are gene segments
//
enum class VJGeneType
{
    IGHV,
    IGHJ,
    IGKV,
    IGKJ,
    IGLV,
    IGLJ,
    TRAV,
    TRAJ,
    TRBV,
    TRBJ,
    TRDV,
    TRDJ,
    TRGV,
    TRGJ;

    val locus: IgTcrLocus = IgTcrLocus.valueOf(name.take(3))
    val vj: VJ = VJ.valueOf(name[3].toString())
}
