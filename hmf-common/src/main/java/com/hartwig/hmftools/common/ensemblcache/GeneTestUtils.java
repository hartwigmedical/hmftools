package com.hartwig.hmftools.common.ensemblcache;

import static java.lang.Math.max;
import static java.lang.Math.min;

import static com.hartwig.hmftools.common.ensemblcache.TranscriptProteinData.BIOTYPE_PROTEIN_CODING;

import java.util.List;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.fusion.GeneAnnotation;
import com.hartwig.hmftools.common.genome.refgenome.RefGenomeVersion;

public class GeneTestUtils
{
    public static GeneAnnotation createGeneAnnotation(int svId, boolean isStart, final String geneName, String stableId, int strand,
            final String chromosome, int position, int orientation)
    {
        String karyotypeBand = "";

        GeneAnnotation gene = new GeneAnnotation(svId, isStart, geneName, stableId, strand, karyotypeBand);
        gene.setPositionalData(chromosome, position, (byte)orientation);

        return gene;
    }

    // Ensembl data types
    public static EnsemblDataCache createGeneDataCache()
    {
        return new EnsemblDataCache("", RefGenomeVersion.HG19);
    }
    public static EnsemblGeneData createEnsemblGeneData(String geneId, String geneName, String chromosome, int strand, int geneStart, int geneEnd)
    {
        return new EnsemblGeneData(geneId, geneName, chromosome, (byte)strand, geneStart, geneEnd,  "");
    }

    public static void addTransExonData(EnsemblDataCache geneTransCache, final String geneId, List<TranscriptData> transDataList)
    {
        geneTransCache.getTranscriptDataMap().put(geneId, transDataList);
    }

    public static void addGeneData(EnsemblDataCache geneTransCache, final String chromosome, List<EnsemblGeneData> geneDataList)
    {
        geneTransCache.getChrGeneDataMap().put(chromosome, geneDataList);
    }

    public static int getCodingBases(final Integer start, final Integer end)
    {
        if(start != null && end != null)
            return (end - start) + 1;
        return 0;
    }

    public static TranscriptData createTransExons(final String geneId, int transId, byte strand,
            int[] exonStarts, int[] exonEndPhases, int exonLength)
    {
        return createTransExons(geneId, transId, strand, exonStarts, exonEndPhases, exonLength, false);
    }

    public static int[] generateExonStarts(int startBase, int exonCount, int exonLength, int intronLength)
    {
        int[] exonStarts = new int[exonCount];

        for(int i = 0; i < exonCount; ++i)
        {
            exonStarts[i] = startBase + i * (exonLength + intronLength);
        }

        return exonStarts;
    }

    public static String generateTransName(int transId) { return String.format("TRAN%04d", transId); }

    public static TranscriptData createTransExons(final String geneId, int transId, byte strand,
            int[] exonStarts, int[] exonEndPhases, int exonLength, boolean isCanonical)
    {
        if(exonStarts.length == 0 || exonStarts.length != exonEndPhases.length)
            return null;

        int exonCount = exonStarts.length;
        int transStart = exonStarts[0];
        int transEnd = exonStarts[exonCount-1] + exonLength;

        Integer codingStart = null;
        Integer codingEnd = null;

        int[] exonPhases = new int[exonCount];

        // work out phases and coding start & end
        for(int i = 0; i < exonCount; ++i)
        {
            int exonStart = exonStarts[i];

            int exonEndPhase = exonEndPhases[i];

            if(strand == 1)
                exonPhases[i] = i > 0 ? exonEndPhases[i-1] : -1;
            else
                exonPhases[i] = i < exonCount - 1 ? exonEndPhases[i+1] : -1;

            if(codingStart == null && ((strand == 1 && exonEndPhase != -1) || (strand == -1 && exonPhases[i] != -1)))
            {
                codingStart = exonStart + exonLength / 2;
            }
            else if(codingStart != null && codingEnd == null
            && ((strand == 1 && exonEndPhase == -1) || (strand == -1 && exonPhases[i] == -1)))
            {
                codingEnd = exonStart + exonLength / 2;
            }
        }

        TranscriptData transData = new TranscriptData(transId, generateTransName(transId), geneId, isCanonical, strand, transStart, transEnd,
                codingStart, codingEnd, BIOTYPE_PROTEIN_CODING);

        List<ExonData> exons = Lists.newArrayList();

        for(int i = 0; i < exonCount; ++i)
        {
            int exonStart = exonStarts[i];
            int exonEnd = exonStarts[i] + exonLength;
            int exonRank = strand == 1 ? i + 1 : exonCount - i;

            exons.add(new ExonData(transId, exonStart, exonEnd, exonRank, exonPhases[i], exonEndPhases[i]));
        }

        transData.setExons(exons);

        return transData;
    }

    public static TranscriptData createTransExons(final String geneId, int transId, byte strand,
            int[] exonStarts, int exonLength, Integer codingStart, Integer codingEnd, boolean isCanonical, final String biotype)
    {
        if(exonStarts.length == 0 || exonLength <= 0)
            return null;

        int exonCount = exonStarts.length;
        int transStart = exonStarts[0];
        int transEnd = exonStarts[exonCount-1] + exonLength;

        final List<ExonData> exons = Lists.newArrayList();

        boolean hasCodingBases = codingStart != null && codingEnd != null;

        // work out phases based on coding start & end
        boolean inCoding = false;
        boolean finishedCoding = false;
        int lastExonEndPhase = -1;

        if(strand == 1)
        {
            for (int i = 0; i < exonCount; ++i)
            {
                int exonStart = exonStarts[i];
                int exonEnd = exonStart + exonLength;
                int exonRank = i + 1;

                int exonPhase = 0;
                int exonStartPhase = -1;
                int exonEndPhase = -1;
                int exonCodingStart = 0;

                if (hasCodingBases && !finishedCoding)
                {
                    if (!inCoding)
                    {
                        if (codingStart <= exonEnd)
                        {
                            inCoding = true;
                            exonPhase = 0;
                            exonStartPhase = codingStart == exonStart ? 0 : -1;
                            exonCodingStart = codingStart;
                        }
                    }
                    else
                    {
                        exonPhase = lastExonEndPhase;
                        exonStartPhase = lastExonEndPhase;
                        exonCodingStart = exonStart;
                    }

                    if (inCoding)
                    {
                        int exonCodingBases = min(exonEnd, codingEnd) - exonCodingStart + 1;
                        exonEndPhase = (exonPhase + exonCodingBases) % 3;
                        lastExonEndPhase = exonEndPhase;

                        if (codingEnd <= exonEnd)
                        {
                            finishedCoding = true;
                            inCoding = false;
                            exonEndPhase = -1;
                        }
                    }
                }

                exons.add(new ExonData(transId, exonStart, exonEnd, exonRank, exonStartPhase, exonEndPhase));
            }
        }
        else
        {
            for (int i = exonCount-1; i >= 0; --i)
            {
                int exonStart = exonStarts[i];
                int exonEnd = exonStart + exonLength;
                int exonRank = exonCount - i;

                int exonPhase = 0;
                int exonStartPhase = -1;
                int exonEndPhase = -1;
                int exonCodingEnd = 0;

                if (hasCodingBases && !finishedCoding)
                {
                    if (!inCoding)
                    {
                        if (codingEnd >= exonStart)
                        {
                            inCoding = true;
                            exonPhase = 0;
                            exonStartPhase = codingStart == exonStart ? 0 : -1;
                            exonCodingEnd = codingEnd;
                        }
                    }
                    else
                    {
                        exonPhase = lastExonEndPhase;
                        exonStartPhase = lastExonEndPhase;
                        exonCodingEnd = exonEnd;
                    }

                    if (inCoding)
                    {
                        int exonCodingBases = exonCodingEnd - max(exonStart, codingStart) + 1;
                        exonEndPhase =  (exonPhase + exonCodingBases) % 3;
                        lastExonEndPhase = exonEndPhase;

                        if (codingStart >= exonStart)
                        {
                            finishedCoding = false;
                            inCoding = false;
                            exonEndPhase = -1;
                        }
                    }
                }

                exons.add(0, new ExonData(transId, exonStart, exonEnd, exonRank, exonStartPhase, exonEndPhase));
            }
        }

        TranscriptData transData = new TranscriptData(transId, generateTransName(transId), geneId, isCanonical, strand, transStart, transEnd,
                codingStart, codingEnd, biotype);

        transData.setExons(exons);

        return transData;

    }

}
