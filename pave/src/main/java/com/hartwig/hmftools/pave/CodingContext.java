package com.hartwig.hmftools.pave;

import static com.hartwig.hmftools.common.gene.CodingBaseData.PHASE_1;
import static com.hartwig.hmftools.common.gene.CodingBaseData.PHASE_2;
import static com.hartwig.hmftools.common.gene.CodingBaseData.PHASE_NONE;
import static com.hartwig.hmftools.common.gene.TranscriptUtils.calcCodingBases;
import static com.hartwig.hmftools.common.gene.TranscriptUtils.getCodingBaseRanges;

import java.util.List;

import com.hartwig.hmftools.common.codon.Codons;
import com.hartwig.hmftools.common.codon.Nucleotides;
import com.hartwig.hmftools.common.gene.CodingBaseData;
import com.hartwig.hmftools.common.gene.TranscriptData;
import com.hartwig.hmftools.common.genome.refgenome.RefGenomeInterface;

// TO-DO - move to hmf-common

public class CodingContext
{
    public int CodingPhase; // phase of first coding base
    public int CodingBase; // index of first coding base for this variant

    public String UpstreamAA;
    public String WildtypeAA;
    public String NovelAA;
    public String DownstreamAA;

    public int BasesToLastExonJunction;

    public CodingContext(
            final int codingPhase, final int codingBase, final String upstreamAA, final String wildtypeAA,
            final String novelAA, final String downstreamAA, final int basesToLastExonJunction)
    {
        CodingPhase = codingPhase;
        CodingBase = codingBase;
        UpstreamAA = upstreamAA;
        WildtypeAA = wildtypeAA;
        NovelAA = novelAA;
        DownstreamAA = downstreamAA;
        BasesToLastExonJunction = basesToLastExonJunction;
    }

    public CodingContext()
    {
        CodingBase = 0;
        CodingPhase = PHASE_NONE;
        UpstreamAA = "";
        WildtypeAA = "";
        NovelAA = "";
        DownstreamAA = "";
        BasesToLastExonJunction = 0;
    }

    public String hgvsCodingChange() { return ""; }
    public String hgvsProteinChange() { return ""; }

    public boolean hasCodingBases() { return !WildtypeAA.isEmpty(); }

    public boolean hasProteinChange() { return !WildtypeAA.equals(NovelAA); }

    public static CodingContext determineContext(
            final String chromosome, int position, final String ref, final String alt,
            final TranscriptData transData, final RefGenomeInterface refGenome)
    {
        // find the start of the codon in which the first base of the variant is in, and likewise the end of the last codon
        int varLength = ref.length();
        boolean posStrand = transData.posStrand();

        int upstreamStartPos = posStrand ? position : position + varLength - 1;
        final CodingBaseData cbData = calcCodingBases(transData, upstreamStartPos);

        // this is 0 if the variant starts on the first base of a codon (phase=1), 1 if it starts on the 2nd, and 2 if on the 3rd/last
        int openCodonBases = getOpenCodonBases(cbData.Phase);

        int codingBaseLen = (int)(Math.ceil((openCodonBases + varLength) / 3.0)) * 3;

        int codonStartPos = posStrand ? upstreamStartPos - openCodonBases : upstreamStartPos + openCodonBases;

        List<int[]> codingRanges = getCodingBaseRanges(transData, codonStartPos, posStrand, codingBaseLen);
        String refCodingBases = refGenome.getBaseString(chromosome, codingRanges);

        String altCodingBases;
        String refAminoAcids;
        String altAminoAcids;

        if(posStrand)
        {
            if(openCodonBases > 0)
            {
                altCodingBases = refCodingBases.substring(0, openCodonBases) + alt
                        + refCodingBases.substring(openCodonBases + varLength);
            }
            else
            {
                altCodingBases = alt + refCodingBases.substring(varLength);
            }

            refAminoAcids = Codons.aminoAcidFromBases(refCodingBases);
            altAminoAcids = Codons.aminoAcidFromBases(altCodingBases);
        }
        else
        {
            if(openCodonBases > 0)
            {
                altCodingBases = refCodingBases.substring(0, codingBaseLen - varLength - openCodonBases) + alt
                        + refCodingBases.substring(codingBaseLen - openCodonBases);
            }
            else
            {
                altCodingBases = refCodingBases.substring(0, codingBaseLen - varLength) + alt;
            }

            refAminoAcids = Codons.aminoAcidFromBases(Nucleotides.reverseStrandBases(refCodingBases));
            altAminoAcids = Codons.aminoAcidFromBases(Nucleotides.reverseStrandBases(altCodingBases));
        }

        int basesToLastExonJunction = 0;
        String upstreamAA = "";
        String downstreamAA = "";

        return new CodingContext(
                cbData.Phase, cbData.CodingBases, upstreamAA, refAminoAcids, altAminoAcids, downstreamAA, basesToLastExonJunction);
    }

    public static int getOpenCodonBases(int phase)
    {
        if(phase == PHASE_1)
            return 0;
        else if(phase == PHASE_2)
            return 1;
        else
            return 2;
    }

}
