package com.hartwig.hmftools.svtools.neo;

import static java.lang.Math.max;
import static java.lang.Math.min;

import static com.hartwig.hmftools.common.fusion.FusionCommon.FS_UPSTREAM;
import static com.hartwig.hmftools.common.fusion.TranscriptCodingType.UTR_3P;
import static com.hartwig.hmftools.common.utils.sv.SvCommonUtils.NEG_ORIENT;
import static com.hartwig.hmftools.common.utils.sv.SvRegion.positionWithin;
import static com.hartwig.hmftools.svtools.neo.AminoAcidConverter.convertDnaCodonToAminoAcid;
import static com.hartwig.hmftools.svtools.neo.AminoAcidConverter.isStopCodon;

import java.util.List;

import com.hartwig.hmftools.common.ensemblcache.ExonData;
import com.hartwig.hmftools.common.ensemblcache.TranscriptData;
import com.hartwig.hmftools.common.fusion.GeneAnnotation;
import com.hartwig.hmftools.common.fusion.Transcript;
import com.hartwig.hmftools.common.fusion.TranscriptRegionType;
import com.hartwig.hmftools.common.genome.refgenome.RefGenomeInterface;

public class NeoUtils
{
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

    public static String getCodingBases(
            final RefGenomeInterface refGenome, final NeoEpitopeData neData, int stream,
            int requiredAminoAcids, int phaseOffset)
    {
        boolean canStartInExon = stream == FS_UPSTREAM || neData.RegionType[FS_UPSTREAM] == TranscriptRegionType.EXONIC;
        int requiredBases = requiredAminoAcids * 3 + phaseOffset;

        return getCodingBases(
                refGenome, neData.TransData[stream], neData.chromosome(stream), neData.positon(stream), neData.orientation(stream),
                requiredBases, canStartInExon);
    }

    public static String getCodingBases(
            final RefGenomeInterface refGenome, final TranscriptData transData,
            final String chromosome, int nePosition, byte neOrientation,
            int requiredBases, boolean canStartInExon)
    {
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


    public static String getBaseStringOld(
            final RefGenomeInterface refGenome, final NeoEpitopeData neData, int stream,
            int requiredAminoAcids, boolean collectAllBases, int phaseOffset)
    {
        final TranscriptData transData = neData.TransData[stream];

        if(transData.CodingStart == null) // TO-DO - get bases to end of transcript anyway for 3' non-coding?
            return "";

        int nePosition = neData.positon(stream);
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

    public static int calcNonMediatedDecayBases(final NeoEpitopeData neData, int stream)
    {
        final List<ExonData> exonDataList = neData.TransData[stream].exons();

        int nePosition = neData.positon(stream);
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
}
