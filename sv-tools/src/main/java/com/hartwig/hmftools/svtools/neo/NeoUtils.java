package com.hartwig.hmftools.svtools.neo;

import static java.lang.Math.max;
import static java.lang.Math.min;

import static com.hartwig.hmftools.svtools.neo.AminoAcidConverter.convertDnaCodonToAminoAcid;
import static com.hartwig.hmftools.svtools.neo.AminoAcidConverter.isStopCodon;
import static com.hartwig.hmftools.svtools.neo.NeoConfig.AMINO_ACID_REF_COUNT;

import java.util.List;

import com.hartwig.hmftools.common.ensemblcache.ExonData;
import com.hartwig.hmftools.common.ensemblcache.TranscriptData;
import com.hartwig.hmftools.common.fusion.GeneAnnotation;
import com.hartwig.hmftools.common.fusion.Transcript;
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

    public static String getBaseString(
            final RefGenomeInterface refGenome, final Transcript transcript, final TranscriptData transData,
            int requiredAminoAcids, boolean collectAllBases, int phaseOffset)
    {
        if(transcript.nonCoding())
            return "";

        final GeneAnnotation gene = transcript.gene();
        int breakPosition = gene.position();

        int requiredBases = requiredAminoAcids * 3 + phaseOffset;

        int codingStart = transcript.CodingStart;
        int codingEnd = transcript.CodingEnd;
        boolean postCoding = transcript.postCoding();

        final List<ExonData> exonDataList = transData.exons();

        String baseString = "";

        if(gene.orientation() == -1)
        {
            // walk forwards through the exons, collecting up the required positions and bases
            for (int i = 0; i < exonDataList.size(); ++i)
            {
                final ExonData exon = exonDataList.get(i);

                if (exon.ExonStart < breakPosition || (!postCoding && exon.ExonStart < codingStart))
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

                baseString += refGenome.getBaseString(gene.chromosome(), posStart, posEnd);

                if (requiredBases <= 0)
                    break;
            }
        }
        else
        {
            for(int i = exonDataList.size() - 1; i >= 0; --i)
            {
                final ExonData exon = exonDataList.get(i);

                if(exon.ExonEnd > breakPosition || (!postCoding && exon.ExonEnd > codingEnd))
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
                baseString = refGenome.getBaseString(gene.chromosome(), posStart, posEnd) + baseString;

                if(requiredBases <= 0)
                    break;
            }

        }

        return baseString;
    }

    public static int calcNonMediatedDecayBases(final GeneAnnotation gene, final TranscriptData transData)
    {
        final List<ExonData> exonDataList = transData.exons();

        int breakPosition = gene.position();

        int exonicBaseCount = 0;

        if(gene.orientation() == -1)
        {
            // walk forwards through the exons, collecting up the required positions and bases
            for (int i = 0; i < exonDataList.size(); ++i)
            {
                final ExonData exon = exonDataList.get(i);

                if(i == exonDataList.size() - 1)
                    break;

                if (breakPosition > exon.ExonEnd)
                    continue;

                exonicBaseCount += exon.ExonEnd - max(breakPosition, exon.ExonStart) + 1;
            }
        }
        else
        {
            for(int i = exonDataList.size() - 1; i >= 0; --i)
            {
                final ExonData exon = exonDataList.get(i);

                if(i == 0)
                    break;

                if(breakPosition < exon.ExonStart)
                    continue;

                exonicBaseCount += min(breakPosition, exon.ExonEnd) - exon.ExonStart + 1;
            }
        }

        return exonicBaseCount;
    }
}
