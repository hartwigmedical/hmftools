package com.hartwig.hmftools.svtools.neo;

import static java.lang.Math.max;
import static java.lang.Math.min;

import static com.hartwig.hmftools.common.fusion.FusionCommon.FS_DOWNSTREAM;
import static com.hartwig.hmftools.common.fusion.FusionCommon.FS_PAIR;
import static com.hartwig.hmftools.common.fusion.FusionCommon.FS_UPSTREAM;
import static com.hartwig.hmftools.common.fusion.FusionCommon.NEG_STRAND;
import static com.hartwig.hmftools.common.fusion.FusionCommon.POS_STRAND;
import static com.hartwig.hmftools.common.fusion.Transcript.CODING_BASES;
import static com.hartwig.hmftools.common.fusion.Transcript.TOTAL_CODING_BASES;
import static com.hartwig.hmftools.common.fusion.TranscriptCodingType.CODING;
import static com.hartwig.hmftools.common.fusion.TranscriptCodingType.NON_CODING;
import static com.hartwig.hmftools.common.fusion.TranscriptCodingType.UTR_3P;
import static com.hartwig.hmftools.common.fusion.TranscriptCodingType.UTR_5P;
import static com.hartwig.hmftools.common.fusion.TranscriptRegionType.EXONIC;
import static com.hartwig.hmftools.common.fusion.TranscriptRegionType.INTRONIC;
import static com.hartwig.hmftools.common.utils.sv.SvCommonUtils.NEG_ORIENT;
import static com.hartwig.hmftools.common.utils.sv.SvRegion.positionWithin;
import static com.hartwig.hmftools.svtools.neo.AminoAcidConverter.convertDnaCodonToAminoAcid;
import static com.hartwig.hmftools.svtools.neo.AminoAcidConverter.isStopCodon;
import static com.hartwig.hmftools.svtools.neo.AminoAcidConverter.reverseStrandBases;

import java.util.List;

import com.hartwig.hmftools.common.ensemblcache.ExonData;
import com.hartwig.hmftools.common.ensemblcache.TranscriptData;
import com.hartwig.hmftools.common.fusion.Transcript;
import com.hartwig.hmftools.common.fusion.TranscriptRegionType;
import com.hartwig.hmftools.common.genome.refgenome.RefGenomeInterface;

public class NeoUtils
{
    public static void setTranscriptContext(
            final NeoEpitope neData, final TranscriptData transData, int position, int stream)
    {
        // determine phasing, coding and region context
        boolean isUpstream = stream == FS_UPSTREAM;

        for(ExonData exon : transData.exons())
        {
            if(position > exon.ExonEnd)
                continue;

            if(position < exon.ExonStart)
            {
                // intronic
                neData.RegionType[stream] = INTRONIC;

                // upstream pos-strand, before next exon then take the exon before's rank
                if(transData.Strand == POS_STRAND && isUpstream)
                    neData.ExonRank[stream] = exon.ExonRank - 1;
                else if(transData.Strand == NEG_STRAND && isUpstream)
                    neData.ExonRank[stream] = exon.ExonRank;
                else if(transData.Strand == POS_STRAND && !isUpstream)
                    neData.ExonRank[stream] = exon.ExonRank;
                else if(transData.Strand == NEG_STRAND && !isUpstream)
                    neData.ExonRank[stream] = exon.ExonRank - 1;

                if(transData.Strand == POS_STRAND)
                    neData.Phases[stream] = exon.ExonPhase;
                else
                    neData.Phases[stream] = exon.ExonPhaseEnd;
            }
            else if(positionWithin(position, exon.ExonStart, exon.ExonEnd))
            {
                neData.RegionType[stream] = EXONIC;
                neData.ExonRank[stream] = exon.ExonRank;
            }

            break;
        }
    }

    public static void setTranscriptCodingData(
            final NeoEpitope neData, final TranscriptData transData, int position, int insSeqLength, int stream)
    {
        if(transData.CodingStart != null)
        {
            if(positionWithin(position, transData.CodingStart, transData.CodingEnd))
                neData.CodingType[stream] = CODING;
            else if(transData.Strand == POS_STRAND && position < transData.CodingStart)
                neData.CodingType[stream] = UTR_5P;
            else if(transData.Strand == NEG_STRAND && position > transData.CodingEnd)
                neData.CodingType[stream] = UTR_5P;
            else
                neData.CodingType[stream] = UTR_3P;
        }
        else
        {
            neData.CodingType[stream] = NON_CODING;
            neData.Phases[stream] = -1;
        }

        if(neData.CodingType[stream] == CODING && neData.RegionType[stream] == EXONIC)
        {
            int[] codingData = Transcript.calcCodingBases(transData.CodingStart, transData.CodingEnd, transData.exons(), position);
            int codingBases = transData.Strand == POS_STRAND ? codingData[CODING_BASES] : codingData[TOTAL_CODING_BASES] - codingData[CODING_BASES] + 1;

            codingBases -= 1;

            // factor in insert sequence for the upstream partner
            codingBases += insSeqLength;

            neData.Phases[stream] = codingBases % 3;
        }
    }

    public static String getCodingBases(
            final RefGenomeInterface refGenome, final TranscriptData transData,
            final String chromosome, int nePosition, byte neOrientation,
            int requiredBases, boolean canStartInExon)
    {
        if(requiredBases <= 0)
            return "";

        int codingStart = transData.CodingStart != null ? transData.CodingStart : transData.TransStart;
        int codingEnd = transData.CodingEnd != null ? transData.CodingEnd : transData.TransEnd;

        final List<ExonData> exonDataList = transData.exons();

        String baseString = "";

        if(neOrientation == NEG_ORIENT)
        {
            // walk forwards through the exons, collecting up the required positions and bases
            for (int i = 0; i < exonDataList.size(); ++i)
            {
                final ExonData exon = exonDataList.get(i);

                if(nePosition > exon.ExonEnd)
                    continue;

                if(codingStart > exon.ExonEnd) // no coding bases reached yet
                    continue;

                if(positionWithin(nePosition, exon.ExonStart, exon.ExonEnd) && !canStartInExon)
                    continue; // will start at the next exon

                if(exon.ExonStart > codingEnd)
                    break; // no more coding bases

                int exonCodingStart = max(codingStart, exon.ExonStart);
                exonCodingStart = max(exonCodingStart, nePosition);
                int exonCodingEnd = min(codingEnd, exon.ExonEnd);
                int exonCodingBases = exonCodingEnd - exonCodingStart + 1;

                int baseStart, baseEnd;

                if(requiredBases >= exonCodingBases)
                {
                    // take them all
                    baseStart = exonCodingStart;
                    baseEnd = exonCodingEnd;
                    requiredBases -= exonCodingBases;
                }
                else
                {
                    baseStart = exonCodingStart;
                    baseEnd = baseStart + requiredBases - 1;
                    requiredBases = 0;
                }

                baseString += refGenome.getBaseString(chromosome, baseStart, baseEnd);

                if (requiredBases <= 0)
                    break;
            }
        }
        else
        {
            for(int i = exonDataList.size() - 1; i >= 0; --i)
            {
                final ExonData exon = exonDataList.get(i);

                if(nePosition < exon.ExonStart)
                    continue;

                if(codingEnd < exon.ExonStart)
                    continue;

                if(positionWithin(nePosition, exon.ExonStart, exon.ExonEnd) && !canStartInExon)
                    continue;

                if(exon.ExonEnd < codingStart)
                    break;

                int exonCodingEnd = min(codingEnd, exon.ExonEnd);
                exonCodingEnd = min(exonCodingEnd, nePosition);
                int exonCodingStart = max(codingStart, exon.ExonStart);
                int exonCodingBases = exonCodingEnd - exonCodingStart + 1;

                int baseStart, baseEnd;

                if(requiredBases >= exonCodingBases)
                {
                    // take them all
                    baseStart = exonCodingStart;
                    baseEnd = exonCodingEnd;
                    requiredBases -= exonCodingBases;
                }
                else
                {
                    baseEnd = exonCodingEnd;
                    baseStart = baseEnd - requiredBases + 1;
                    requiredBases = 0;
                }

                baseString = refGenome.getBaseString(chromosome, baseStart, baseEnd) + baseString;

                if (requiredBases <= 0)
                    break;
            }
        }

        return baseString;
    }

    public static void adjustCodingBasesForStrand(final NeoEpitope neData)
    {
        // upstream strand 1, bases will be retrieved from left to right (lower to higher), no need for any conversion
        // downstream strand 1, bases will be retrieved from left to right (lower to higher), no need for any conversion
        // upstream strand -1, bases will be retrieved from left to right (lower to higher), need to reverse and convert
        // downstream strand -1, bases will be retrieved from left to right (lower to higher), need to reverse and convert

        for(int fs = FS_UPSTREAM; fs <= FS_DOWNSTREAM; ++fs)
        {
            if(neData.strand(fs) == NEG_STRAND)
                neData.CodingBases[fs] = reverseStrandBases(neData.CodingBases[fs]);
        }
    }

    public static String getBaseStringOld(
            final RefGenomeInterface refGenome, final NeoEpitope neData, int stream,
            int requiredAminoAcids, boolean collectAllBases, int phaseOffset)
    {
        final TranscriptData transData = neData.TransData[stream];

        if(transData.CodingStart == null) // TO-DO - get bases to end of transcript anyway for 3' non-coding?
            return "";

        int nePosition = neData.position(stream);
        final String chromosome = neData.chromosome(stream);
        int requiredBases = requiredAminoAcids * 3 + phaseOffset;

        int codingStart = transData.CodingStart;
        int codingEnd = transData.CodingEnd;

        boolean postCoding = neData.CodingType[stream] == UTR_3P;

        final List<ExonData> exonDataList = transData.exons();

        String baseString = "";

        if(neData.orientation(stream) == NEG_ORIENT)
        {
            // walk forwards through the exons, collecting up the required positions and bases
            for (int i = 0; i < exonDataList.size(); ++i)
            {
                final ExonData exon = exonDataList.get(i);

                if (exon.ExonStart < nePosition || (!postCoding && exon.ExonStart < codingStart))
                    continue;

                int posStart, posEnd;

                if (collectAllBases || requiredBases > exon.ExonEnd - exon.ExonStart + 1)
                {
                    posStart = exon.ExonStart;
                    posEnd = exon.ExonEnd;
                    requiredBases -= (exon.ExonEnd - exon.ExonStart + 1);
                }
                else
                {
                    posStart = exon.ExonStart;
                    posEnd = exon.ExonStart + requiredBases - 1;
                    requiredBases = 0;
                }

                // stop at end of coding region unless the breakend started past it
                if(!postCoding && posEnd > codingEnd)
                {
                    posEnd = codingEnd;
                    requiredBases = 0;
                }

                if(posEnd < posStart)
                    continue;

                baseString += refGenome.getBaseString(chromosome, posStart, posEnd);

                if (requiredBases <= 0)
                    break;
            }
        }
        else
        {
            for(int i = exonDataList.size() - 1; i >= 0; --i)
            {
                final ExonData exon = exonDataList.get(i);

                if(exon.ExonEnd > nePosition || (!postCoding && exon.ExonEnd > codingEnd))
                    continue;

                int posStart, posEnd;

                if(collectAllBases || requiredBases > exon.ExonEnd - exon.ExonStart + 1)
                {
                    posEnd = exon.ExonEnd;
                    posStart = exon.ExonStart;
                    requiredBases -= (exon.ExonEnd - exon.ExonStart + 1);
                }
                else
                {
                    posEnd = exon.ExonEnd;
                    posStart = exon.ExonEnd - requiredBases + 1;
                    requiredBases = 0;
                }

                if(!postCoding && posStart < codingStart)
                {
                    posStart = codingStart;
                    requiredBases = 0;
                }

                if(posEnd < posStart)
                    continue;

                // add in reverse since walking backwards through the exons
                baseString = refGenome.getBaseString(chromosome, posStart, posEnd) + baseString;

                if(requiredBases <= 0)
                    break;
            }
        }

        return baseString;
    }

    public static int calcNonMediatedDecayBases(final NeoEpitope neData, int stream)
    {
        final List<ExonData> exonDataList = neData.TransData[stream].exons();

        int nePosition = neData.position(stream);
        int exonicBaseCount = 0;

        if(neData.orientation(stream) == -1)
        {
            // walk forwards through the exons, collecting up the required positions and bases
            for (int i = 0; i < exonDataList.size(); ++i)
            {
                final ExonData exon = exonDataList.get(i);

                if(i == exonDataList.size() - 1)
                    break;

                if (nePosition > exon.ExonEnd)
                    continue;

                exonicBaseCount += exon.ExonEnd - max(nePosition, exon.ExonStart) + 1;
            }
        }
        else
        {
            for(int i = exonDataList.size() - 1; i >= 0; --i)
            {
                final ExonData exon = exonDataList.get(i);

                if(i == 0)
                    break;

                if(nePosition < exon.ExonStart)
                    continue;

                exonicBaseCount += min(nePosition, exon.ExonEnd) - exon.ExonStart + 1;
            }
        }

        return exonicBaseCount;
    }

    public static String checkTrimBases(final String bases)
    {
        if(bases.length() < 50)
            return bases;

        return bases.substring(0, 50) + "...";
    }

    public static String getAminoAcids(final String baseString, boolean checkStopCodon)
    {
        if(baseString.length() < 3)
            return "";

        String aminoAcidStr = "";
        int index = 0;
        while(index <= baseString.length() - 3)
        {
            String codonBases = baseString.substring(index, index + 3);

            if(isStopCodon(codonBases) && checkStopCodon)
                break;

            String aminoAcid = convertDnaCodonToAminoAcid(codonBases);

            aminoAcidStr += aminoAcid;
            index += 3;
        }

        return aminoAcidStr;
    }

}
