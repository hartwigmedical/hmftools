package com.hartwig.hmftools.pave;

import static java.lang.Math.max;
import static java.lang.Math.min;

import static com.hartwig.hmftools.common.gene.CodingBaseData.PHASE_1;
import static com.hartwig.hmftools.common.gene.CodingBaseData.PHASE_2;
import static com.hartwig.hmftools.common.gene.TranscriptUtils.calcCodingBases;
import static com.hartwig.hmftools.common.gene.TranscriptUtils.getCodingBaseRanges;
import static com.hartwig.hmftools.common.utils.sv.StartEndIterator.SE_END;
import static com.hartwig.hmftools.common.utils.sv.StartEndIterator.SE_START;

import java.util.List;

import com.hartwig.hmftools.common.codon.Codons;
import com.hartwig.hmftools.common.codon.Nucleotides;
import com.hartwig.hmftools.common.gene.CodingBaseData;
import com.hartwig.hmftools.common.gene.ExonData;
import com.hartwig.hmftools.common.gene.TranscriptData;
import com.hartwig.hmftools.common.genome.refgenome.RefGenomeInterface;

public class ProteinContext
{
    public String UpstreamAA;
    public String WildtypeAA;
    public String NovelAA;
    public String DownstreamAA;

    public ProteinContext()
    {
        UpstreamAA = "";
        WildtypeAA = "";
        NovelAA = "";
        DownstreamAA = "";
    }

    public boolean hasCodingBases() { return !WildtypeAA.isEmpty(); }

    public boolean hasProteinChange() { return !WildtypeAA.equals(NovelAA); }

    public static ProteinContext determineContext(
            final VariantData variant, final CodingContext codingContext, final TranscriptData transData, final RefGenomeInterface refGenome)
    {
        // find the start of the codon in which the first base of the variant is in, and likewise the end of the last codon
        int varLength = variant.Ref.length();
        boolean posStrand = transData.posStrand();

        // int upstreamStartPos = posStrand ? variant.Position : variant.Position + varLength - 1;
        // final CodingBaseData cbData = calcCodingBases(transData, upstreamStartPos);

        // this is 0 if the variant starts on the first base of a codon (phase=1), 1 if it starts on the 2nd, and 2 if on the 3rd/last
        int upstreamStartPos = posStrand ? variant.Position : variant.EndPosition;
        int openCodonBases = getOpenCodonBases(codingContext.CodingPhase);

        int codingBaseLen = (int)(Math.ceil((openCodonBases + varLength) / 3.0)) * 3;

        int codonStartPos = posStrand ? upstreamStartPos - openCodonBases : upstreamStartPos + openCodonBases;

        int stream = posStrand ? SE_START : SE_END;
        ExonData exon = transData.exons().stream().filter(x -> x.Rank == codingContext.ExonRank[stream]).findFirst().orElse(null);

        int codingPosStart;
        int codingPosEnd;

        if(posStrand)
        {
            codingPosStart = max(codonStartPos, exon.Start);
            codingPosEnd = min(codonStartPos + codingBaseLen - 1, transData.CodingEnd);
        }
        else
        {
            codingPosEnd = min(codonStartPos, exon.End);
            codingPosStart = max(codonStartPos - codingBaseLen + 1, transData.CodingStart);
        }

        // List<int[]> codingRanges = getCodingBaseRanges(transData, codonStartPos, posStrand, codingBaseLen);
        String refCodingBases = refGenome.getBaseString(variant.Chromosome, codingPosStart, codingPosEnd);

        ProteinContext pc = new ProteinContext();

        String alt = variant.Alt;
        String altCodingBases;

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

            pc.WildtypeAA = Codons.aminoAcidFromBases(refCodingBases);
            pc.NovelAA = Codons.aminoAcidFromBases(altCodingBases);
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

            pc.WildtypeAA = Codons.aminoAcidFromBases(Nucleotides.reverseStrandBases(refCodingBases));
            pc.NovelAA = Codons.aminoAcidFromBases(Nucleotides.reverseStrandBases(altCodingBases));
        }

        // pc.UpstreamAA = ;
        // pc.DownstreamAA = ;

        return pc;
    }

    private static int getOpenCodonBases(int phase)
    {
        if(phase == PHASE_1)
            return 0;
        else if(phase == PHASE_2)
            return 1;
        else
            return 2;
    }

}
