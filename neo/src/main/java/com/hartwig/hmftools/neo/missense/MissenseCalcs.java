package com.hartwig.hmftools.neo.missense;

import static java.lang.Math.max;
import static java.lang.Math.min;

import static com.hartwig.hmftools.common.codon.Codons.STOP_AMINO_ACID;
import static com.hartwig.hmftools.common.codon.Codons.isStopCodon;
import static com.hartwig.hmftools.common.codon.Nucleotides.reverseComplementBases;
import static com.hartwig.hmftools.common.codon.Nucleotides.swapDnaBase;
import static com.hartwig.hmftools.common.fusion.FusionCommon.POS_STRAND;
import static com.hartwig.hmftools.common.gene.TranscriptUtils.calcCodingStartPositionAdjustment;
import static com.hartwig.hmftools.common.utils.sv.SvCommonUtils.POS_ORIENT;
import static com.hartwig.hmftools.neo.NeoCommon.NE_LOGGER;

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
                    exonCodingStart += calcCodingStartPositionAdjustment(transData, exon);
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
                    exonCodingEnd += calcCodingStartPositionAdjustment(transData, exon);
                }

                codingBases = mRefGenome.getBaseString(geneData.Chromosome, exonCodingStart, exonCodingEnd) + codingBases;

                for(int pos = exonCodingEnd; pos >= exonCodingStart; --pos)
                {
                    codingBasePositions.add(pos);
                }
            }

            produceMissensePeptides(geneData, transData, reverseComplementBases(codingBases), codingBasePositions);
        }
    }

    private void produceMissensePeptides(
            final GeneData geneData, final TranscriptData transData, final String codingBases, final List<Integer> codingBasePositions)
    {
        int codonCount = codingBases.length() / 3;

        MissensePeptide peptideData = new MissensePeptide(geneData.GeneId, geneData.GeneName, transData.TransName);

        // create peptides for the specified length range and with flanks
        for(int codonIndex = 0; codonIndex < codonCount; ++codonIndex)
        {
            peptideData.CodonIndex = codonIndex + 1; // 1 based
            int codonStartBaseIndex = codonIndex * 3;

            // cycle through each missense variant
            String refCodon = codingBases.substring(codonStartBaseIndex, codonStartBaseIndex + 3);

            String codonRefAminoAcid = isStopCodon(refCodon) ? String.valueOf(STOP_AMINO_ACID) : AminoAcids.findAminoAcidForCodon(refCodon);

            if(codonRefAminoAcid == null)
                return;

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

                    peptideData.Context = geneData.Strand == POS_ORIENT ? refCodon : reverseComplementBases(refCodon);
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
        for(int pepLen = mConfig.PeptideLengthMin; pepLen <= mConfig.PeptideLengthMax; ++pepLen)
        {
            generatePeptides(basePeptideData, codingBases, codonIndex, altAminoAcid, pepLen);
        }
    }

    private void generatePeptides(
            final MissensePeptide basePeptideData, final String codingBases, int targetCodonIndex, final String altAminoAcid, int pepLen)
    {
        // if say length is 8 and the current codon index is 15 with flank of 3, then the peptide + flanks needs to cover
        // the codon range 8 -> 15 to 15 -> 22, which flanks becomes 5 -> 15 up to 15 -> 25
        // where enough AA / codons are missing flanks, take whatever is available
        int peptideCodonStartIndex = max(targetCodonIndex - pepLen + 1, 0);

        NE_LOGGER.trace("gene({}) pepLen({}) codonIndex({}) codonRangeStart({})",
                basePeptideData.GeneName, pepLen, targetCodonIndex, peptideCodonStartIndex);

        int codonCount = codingBases.length() / 3;

        // test peptides of the specified length

        // for(int pepStartCodonIndex = peptideCodonStartIndex; pepStartCodonIndex <= codonIndex; ++pepStartCodonIndex)
        for(int codonIndex = peptideCodonStartIndex; codonIndex <= targetCodonIndex; ++codonIndex)
        {
            int peptideEndCodingBaseIndex = (codonIndex + pepLen - 1) * 3;
            if(peptideEndCodingBaseIndex >= codingBases.length())
                break;

            String peptide = "";

            for(int pepIndex = 0; pepIndex < pepLen; ++pepIndex)
            {
                int codingBaseCodonIndex = codonIndex + pepIndex;

                if(codingBaseCodonIndex == targetCodonIndex)
                {
                    peptide += altAminoAcid;
                }
                else
                {
                    int codonBaseIndex = (codonIndex + pepIndex) * 3;

                    String aminoAcid = null;
                    if(codonBaseIndex + 3 <= codingBases.length())
                    {
                        String codon = codingBases.substring(codonBaseIndex, codonBaseIndex + 3);
                        aminoAcid = AminoAcids.findAminoAcidForCodon(codon);
                    }
                    else
                    {
                        NE_LOGGER.warn("gene({}) codon request({}-{}) overruns coding bases length({})",
                                basePeptideData.GeneName, codonBaseIndex, codonBaseIndex + 3, codingBases.length());
                    }

                    if(aminoAcid == null)
                    {
                        peptide = "";
                        break;
                    }

                    peptide += aminoAcid;
                }
            }

            if(peptide.isEmpty())
                continue;

            MissensePeptide peptideData = new MissensePeptide(basePeptideData.GeneId, basePeptideData.GeneName, basePeptideData.TransName);
            peptideData.Position = basePeptideData.Position;
            peptideData.CodonIndex = basePeptideData.CodonIndex;
            peptideData.Context = basePeptideData.Context;
            peptideData.RefBase = basePeptideData.RefBase;
            peptideData.AltBase = basePeptideData.AltBase;
            peptideData.Peptide = peptide;
            peptideData.PeptideStartIndex = codonIndex + 1; // also 1-based

            if(!mConfig.KeepDuplicates)
            {
                if(mUniqiuePeptides.contains(peptideData.Peptide))
                    continue;

                mUniqiuePeptides.add(peptideData.Peptide);
            }

            if(mConfig.FlankLength > 0)
            {
                String upFlank = "";

                int flankUpStartIndex = max(codonIndex - mConfig.FlankLength, 0);
                for(int flankIndex = flankUpStartIndex; flankIndex < codonIndex; ++flankIndex)
                {
                    int codonBaseIndex = flankIndex * 3;
                    String codon = codingBases.substring(codonBaseIndex, codonBaseIndex + 3);
                    String aminoAcid = AminoAcids.findAminoAcidForCodon(codon);

                    if(aminoAcid == null)
                    {
                        upFlank = "";
                        break;
                    }

                    upFlank += aminoAcid;
                }

                peptideData.UpFlank = upFlank;

                String downFlank = "";

                int flankDownStartIndex = codonIndex + pepLen;
                int flankDownEndIndex = min(flankDownStartIndex + mConfig.FlankLength, codonCount);
                for(int flankIndex = flankDownStartIndex; flankIndex < flankDownEndIndex; ++flankIndex)
                {
                    int codonBaseIndex = flankIndex * 3;
                    String codon = codingBases.substring(codonBaseIndex, codonBaseIndex + 3);
                    String aminoAcid = AminoAcids.findAminoAcidForCodon(codon);

                    if(aminoAcid == null)
                        break;

                    downFlank += aminoAcid;
                }

                peptideData.DownFlank = downFlank;
            }

            mPeptideData.add(peptideData);
        }
    }
}
