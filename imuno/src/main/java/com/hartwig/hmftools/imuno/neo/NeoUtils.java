package com.hartwig.hmftools.imuno.neo;

import static java.lang.Math.max;
import static java.lang.Math.min;

import static com.hartwig.hmftools.common.fusion.FusionCommon.FS_DOWN;
import static com.hartwig.hmftools.common.fusion.FusionCommon.FS_UP;
import static com.hartwig.hmftools.common.fusion.FusionCommon.NEG_STRAND;
import static com.hartwig.hmftools.common.fusion.FusionCommon.POS_STRAND;
import static com.hartwig.hmftools.common.fusion.TranscriptCodingType.CODING;
import static com.hartwig.hmftools.common.fusion.TranscriptCodingType.NON_CODING;
import static com.hartwig.hmftools.common.fusion.TranscriptCodingType.UTR_3P;
import static com.hartwig.hmftools.common.fusion.TranscriptCodingType.UTR_5P;
import static com.hartwig.hmftools.common.fusion.TranscriptRegionType.EXONIC;
import static com.hartwig.hmftools.common.fusion.TranscriptRegionType.INTRONIC;
import static com.hartwig.hmftools.common.fusion.TranscriptUtils.calcCodingBases;
import static com.hartwig.hmftools.common.fusion.TranscriptUtils.codingBasesToPhase;
import static com.hartwig.hmftools.common.fusion.TranscriptUtils.tickPhaseForward;
import static com.hartwig.hmftools.common.utils.sv.SvCommonUtils.NEG_ORIENT;
import static com.hartwig.hmftools.common.utils.sv.SvRegion.positionWithin;
import static com.hartwig.hmftools.imuno.neo.AminoAcidConverter.convertDnaCodonToAminoAcid;
import static com.hartwig.hmftools.imuno.neo.AminoAcidConverter.isStopCodon;
import static com.hartwig.hmftools.imuno.neo.AminoAcidConverter.reverseStrandBases;

import java.util.List;

import com.hartwig.hmftools.common.ensemblcache.ExonData;
import com.hartwig.hmftools.common.ensemblcache.TranscriptData;
import com.hartwig.hmftools.common.fusion.CodingBaseData;
import com.hartwig.hmftools.common.genome.refgenome.RefGenomeInterface;

public class NeoUtils
{
    public static void setTranscriptContext(
            final NeoEpitope neData, final TranscriptData transData, int position, int stream)
    {
        // determine phasing, coding and region context
        boolean isUpstream = stream == FS_UP;

        for(ExonData exon : transData.exons())
        {
            if(position > exon.End)
                continue;

            if(position < exon.Start)
            {
                // intronic
                neData.RegionType[stream] = INTRONIC;

                // upstream pos-strand, before next exon then take the exon before's rank
                if(transData.Strand == POS_STRAND && isUpstream)
                    neData.ExonRank[stream] = exon.Rank - 1;
                else if(transData.Strand == NEG_STRAND && isUpstream)
                    neData.ExonRank[stream] = exon.Rank;
                else if(transData.Strand == POS_STRAND && !isUpstream)
                    neData.ExonRank[stream] = exon.Rank;
                else if(transData.Strand == NEG_STRAND && !isUpstream)
                    neData.ExonRank[stream] = exon.Rank - 1;

                if(transData.Strand == POS_STRAND)
                    neData.Phases[stream] = exon.PhaseStart;
                else
                    neData.Phases[stream] = exon.PhaseEnd;
            }
            else if(positionWithin(position, exon.Start, exon.End))
            {
                neData.RegionType[stream] = EXONIC;
                neData.ExonRank[stream] = exon.Rank;
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
            final CodingBaseData cbData = calcCodingBases(transData, position);
            neData.Phases[stream] = tickPhaseForward(cbData.Phase, insSeqLength);
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

                if(nePosition > exon.End)
                    continue;

                if(codingStart > exon.End) // no coding bases reached yet
                    continue;

                if(positionWithin(nePosition, exon.Start, exon.End) && !canStartInExon)
                    continue; // will start at the next exon

                if(exon.Start > codingEnd)
                    break; // no more coding bases

                int exonCodingStart = max(codingStart, exon.Start);
                exonCodingStart = max(exonCodingStart, nePosition);
                int exonCodingEnd = min(codingEnd, exon.End);
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

                if(nePosition < exon.Start)
                    continue;

                if(codingEnd < exon.Start)
                    continue;

                if(positionWithin(nePosition, exon.Start, exon.End) && !canStartInExon)
                    continue;

                if(exon.End < codingStart)
                    break;

                int exonCodingEnd = min(codingEnd, exon.End);
                exonCodingEnd = min(exonCodingEnd, nePosition);
                int exonCodingStart = max(codingStart, exon.Start);
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

        for(int fs = FS_UP; fs <= FS_DOWN; ++fs)
        {
            if(neData.strand(fs) == NEG_STRAND)
                neData.CodingBases[fs] = reverseStrandBases(neData.CodingBases[fs]);
        }
    }

    public static int calcNonMediatedDecayBases(final NeoEpitope neData, int stream)
    {
        final TranscriptData transData = neData.TransData[stream];
        final List<ExonData> exonDataList = transData.exons();

        int exonicBaseCount = 0;

        if(neData.orientation(stream) == NEG_STRAND)
        {
            int refPosition = transData.CodingEnd != null ? transData.CodingEnd : neData.position(stream);

            for (int i = 0; i < exonDataList.size(); ++i)
            {
                final ExonData exon = exonDataList.get(i);

                if(i == exonDataList.size() - 1)
                    break;

                if (refPosition > exon.End)
                    continue;

                exonicBaseCount += exon.End - max(refPosition, exon.Start) + 1;
            }
        }
        else
        {
            int refPosition = transData.CodingStart != null ? transData.CodingStart : neData.position(stream);

            for(int i = exonDataList.size() - 1; i >= 0; --i)
            {
                final ExonData exon = exonDataList.get(i);

                if(i == 0)
                    break;

                if(refPosition < exon.Start)
                    continue;

                exonicBaseCount += min(refPosition, exon.End) - exon.Start + 1;
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
