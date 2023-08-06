package com.hartwig.hmftools.neo.missense;

import static java.lang.Math.max;
import static java.lang.Math.min;

import static com.hartwig.hmftools.common.codon.Nucleotides.reverseStrandBases;
import static com.hartwig.hmftools.common.codon.Nucleotides.swapDnaBase;
import static com.hartwig.hmftools.common.fusion.FusionCommon.POS_STRAND;
import static com.hartwig.hmftools.common.gene.CodingBaseData.PHASE_0;
import static com.hartwig.hmftools.common.gene.CodingBaseData.PHASE_2;
import static com.hartwig.hmftools.common.gene.CodingBaseData.PHASE_NONE;
import static com.hartwig.hmftools.common.utils.sv.SvCommonUtils.POS_ORIENT;
import static com.hartwig.hmftools.neo.NeoCommon.NE_LOGGER;
import static com.hartwig.hmftools.neo.bind.BindScorer.INVALID_CALC;

import java.util.List;
import java.util.Set;

import com.google.common.collect.Lists;
import com.google.common.collect.Sets;
import com.hartwig.hmftools.common.codon.AminoAcids;
import com.hartwig.hmftools.common.codon.Nucleotides;
import com.hartwig.hmftools.common.gene.ExonData;
import com.hartwig.hmftools.common.gene.GeneData;
import com.hartwig.hmftools.common.gene.TranscriptData;
import com.hartwig.hmftools.common.genome.refgenome.RefGenomeInterface;

public class MissenseCalcs
{
    private final MissenseConfig mConfig;
    private final RefGenomeInterface mRefGenome;

    private final Set<String> mUniqiuePeptides;
    private final List<MissensePeptide> mPeptideData;
    private final List<MissenseVariant> mMissenseVariants;

    public MissenseCalcs(final MissenseConfig config, final RefGenomeInterface refGenome)
    {
        mConfig = config;
        mRefGenome = refGenome;

        mUniqiuePeptides = Sets.newHashSet();
        mPeptideData = Lists.newArrayList();
        mMissenseVariants = Lists.newArrayList();
    }

    public List<MissensePeptide> peptideData() { return mPeptideData; }

    public void clear()
    {
        mUniqiuePeptides.clear();
        mPeptideData.clear();
        mMissenseVariants.clear();
    }

    public void processTranscript(final GeneData geneData, final TranscriptData transData)
    {
        if(transData.CodingStart == null)
            return;

        boolean inCoding = false;

        List<Integer> codingBasePositions = Lists.newArrayList();

        if(transData.Strand == POS_STRAND)
        {
            StringBuilder codingBases = new StringBuilder();

            for(int i = 0; i < transData.exons().size(); ++i)
            {
                ExonData exon = transData.exons().get(i);

                if(exon.End < transData.CodingStart)
                    continue;

                if(exon.Start > transData.CodingEnd)
                    break;

                int exonCodingStart = max(transData.CodingStart, exon.Start);
                int exonCodingEnd = min(transData.CodingEnd, exon.End);

                if(!inCoding)
                {
                    inCoding = true;

                    if(exon.Start == transData.CodingStart)
                    {
                        int startPhase = exon.PhaseStart == PHASE_NONE ? PHASE_0 : exon.PhaseStart;

                        if(startPhase == PHASE_2)
                            exonCodingStart += 2;
                        else if(startPhase == PHASE_0)
                            exonCodingStart += 1;
                    }
                }

                for(int pos = exonCodingStart; pos <= exonCodingEnd; ++pos)
                {
                    codingBasePositions.add(pos);
                }

                codingBases.append(mRefGenome.getBaseString(geneData.Chromosome, exonCodingStart, exonCodingEnd));
            }

            produceMissensePeptides(geneData, transData, codingBases.toString(), codingBasePositions);
        }
        else
        {
            String codingBases = "";
            for(int i = transData.exons().size() - 1; i >= 0; --i)
            {
                ExonData exon = transData.exons().get(i);

                if(exon.Start > transData.CodingEnd)
                    continue;

                if(exon.End < transData.CodingStart)
                    break;

                int exonCodingStart = max(transData.CodingStart, exon.Start);
                int exonCodingEnd = min(transData.CodingEnd, exon.End);

                if(!inCoding)
                {
                    inCoding = true;

                    if(exon.End == transData.CodingEnd)
                    {
                        int startPhase = exon.PhaseStart == PHASE_NONE ? PHASE_0 : exon.PhaseStart;

                        if(startPhase == PHASE_2)
                            exonCodingStart -= 2;
                        else if(startPhase == PHASE_0)
                            exonCodingStart -= 1;
                    }
                }

                codingBases = mRefGenome.getBaseString(geneData.Chromosome, exonCodingStart, exonCodingEnd) + codingBases;

                for(int pos = exonCodingEnd; pos >= exonCodingStart; --pos)
                {
                    codingBasePositions.add(pos);
                }
            }

            produceMissensePeptides(geneData, transData, reverseStrandBases(codingBases), codingBasePositions);
        }
    }

    private void produceMissensePeptides(
            final GeneData geneData, final TranscriptData transData, final String codingBases, final List<Integer> codingBasePositions)
    {
        int codonCount = codingBases.length() / 3;

        MissensePeptide peptideData = new MissensePeptide(geneData.GeneId, geneData.GeneName, transData.TransName);

        // create peptides for the specified length range and with flanks
        for(int codonIndex = mConfig.FlankLength; codonIndex < codonCount - mConfig.FlankLength; ++codonIndex)
        {
            peptideData.CodonIndex = codonIndex + 1; // 1 based
            int codonStartBaseIndex = codonIndex * 3;

            // cycle through each missense variant
            String refCodon = codingBases.substring(codonStartBaseIndex, codonStartBaseIndex + 3);
            String codonRefAminoAcid = AminoAcids.findAminoAcidForCodon(refCodon);

            for(int codonBaseIndex = 0; codonBaseIndex <= 2; ++codonBaseIndex)
            {
                char codonRefBase = refCodon.charAt(codonBaseIndex);
                peptideData.Position = codingBasePositions.get(codonStartBaseIndex + codonBaseIndex);

                for(char dnaBase : Nucleotides.DNA_BASES)
                {
                    if(dnaBase == codonRefBase)
                        continue;

                    String altCodon = "";

                    for(int i = 0; i <= 2; ++i)
                    {
                        if(i == codonBaseIndex)
                            altCodon += dnaBase;
                        else
                            altCodon += codingBases.charAt(codonIndex * 3 + i);
                    }

                    peptideData.Context = geneData.Strand == POS_ORIENT ? refCodon : reverseStrandBases(refCodon);
                    peptideData.RefBase = geneData.Strand == POS_ORIENT ? codonRefBase : swapDnaBase(codonRefBase);
                    peptideData.AltBase = geneData.Strand == POS_ORIENT ? dnaBase : swapDnaBase(dnaBase);

                    String codonAltAminoAcid = AminoAcids.findAminoAcidForCodon(altCodon);

                    if(codonAltAminoAcid == null || codonRefAminoAcid.equals(codonAltAminoAcid))
                        continue;

                    generatePeptides(peptideData, codingBases, codonIndex, codonAltAminoAcid);
                }
            }
        }
    }

    private void generatePeptides(
            final MissensePeptide basePeptideData, final String codingBases, int codonIndex, final String altAminoAcid)
    {
        for(int pepLen = mConfig.PeptideLengths[0]; pepLen <= mConfig.PeptideLengths[1]; ++pepLen)
        {
            // if say length is 8 and the current codon index is 15 with flank of 3, then the peptide + flanks needs to cover
            // the codon range 8 -> 15 to 15 -> 22, which flanks becomes 5 -> 15 up to 15 -> 25
            int peptideCodonStartIndex = max(0, codonIndex - pepLen - mConfig.FlankLength + 1);

            NE_LOGGER.trace("gene({}) pepLen({}) codonIndex({}) codonRangeStart({})",
                    basePeptideData.GeneName, pepLen, codonIndex, peptideCodonStartIndex);

            for(int pepStartCodonIndex = peptideCodonStartIndex; pepStartCodonIndex <= codonIndex - mConfig.FlankLength; ++pepStartCodonIndex)
            {
                int pepPlusFlanksLen = pepLen + mConfig.FlankLength * 2;

                int peptideEndCodingBaseIndex = (pepStartCodonIndex + pepPlusFlanksLen - 1) * 3;
                if(peptideEndCodingBaseIndex >= codingBases.length())
                    break;

                String peptidePlusFlanks = "";

                for(int pepIndex = 0; pepIndex < pepPlusFlanksLen; ++pepIndex)
                {
                    int codingBaseCodonIndex = pepStartCodonIndex + pepIndex;

                    if(codingBaseCodonIndex == codonIndex)
                    {
                        peptidePlusFlanks += altAminoAcid;
                    }
                    else
                    {
                        int pepCodonBaseIndex = (pepStartCodonIndex + pepIndex) * 3;
                        String codon = codingBases.substring(pepCodonBaseIndex, pepCodonBaseIndex + 3);
                        String aminoAcid = AminoAcids.findAminoAcidForCodon(codon);;

                        if(aminoAcid == null)
                        {
                            peptidePlusFlanks = "";
                            break;
                        }

                        peptidePlusFlanks += aminoAcid;
                    }
                }

                if(peptidePlusFlanks.isEmpty() || peptidePlusFlanks.length() < mConfig.FlankLength * 2 + 1)
                    continue;

                MissensePeptide peptideData = new MissensePeptide(basePeptideData.GeneId, basePeptideData.GeneName, basePeptideData.TransName);
                peptideData.Position = basePeptideData.Position;
                peptideData.CodonIndex = basePeptideData.CodonIndex;
                peptideData.Context = basePeptideData.Context;
                peptideData.RefBase = basePeptideData.RefBase;
                peptideData.AltBase = basePeptideData.AltBase;

                peptideData.Peptide = peptidePlusFlanks.substring(mConfig.FlankLength, peptidePlusFlanks.length() - mConfig.FlankLength);

                if(!mConfig.KeepDuplicates)
                {
                    if(mUniqiuePeptides.contains(peptideData.Peptide))
                        continue;

                    mUniqiuePeptides.add(peptideData.Peptide);
                }

                peptideData.UpFlank = peptidePlusFlanks.substring(0, mConfig.FlankLength);
                peptideData.DownFlank = peptidePlusFlanks.substring(peptidePlusFlanks.length() - mConfig.FlankLength);

                mPeptideData.add(peptideData);
            }
        }
    }

    public boolean passesRankThreshold(double value)
    {
        return value != INVALID_CALC && (mConfig.LikelihoodCutoff == 0 || mConfig.LikelihoodCutoff > 0 && value <= mConfig.LikelihoodCutoff);
    }
}
