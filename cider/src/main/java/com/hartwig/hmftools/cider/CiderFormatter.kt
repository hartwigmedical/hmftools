package com.hartwig.hmftools.cider

import com.hartwig.hmftools.common.codon.Codons

// object to just help us format the VDJ data
// the way we want
object CiderFormatter
{
    // we want to replace X with _, easier to see in file
    fun aminoAcidFromBases(dna: String): String
    {
        return Codons.aminoAcidFromBases(dna).replace(Codons.STOP_AMINO_ACID, '_')
    }

    // add suffix for out of frame
    fun cdr3AminoAcid(vdj: VDJSequence): String
    {
        val suffix = if (vdj.isInFrame) "" else "fs"
        return aminoAcidFromBases(vdj.cdr3Sequence) + suffix
    }

    fun vAnchorAA(vdj: VDJSequence): String
    {
        val vAnchorSeq = vdj.vAnchorSequence

        // codon align
        return aminoAcidFromBases(vAnchorSeq.drop(vAnchorSeq.length % 3))
    }

    fun jAnchorAA(vdj: VDJSequence): String
    {
        return aminoAcidFromBases(vdj.jAnchorSequence)
    }
}