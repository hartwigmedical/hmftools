package com.hartwig.hmftools.linx.gene;

import java.util.List;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.variant.structural.annotation.EnsemblGeneData;
import com.hartwig.hmftools.common.variant.structural.annotation.ExonData;
import com.hartwig.hmftools.common.variant.structural.annotation.GeneAnnotation;
import com.hartwig.hmftools.common.variant.structural.annotation.TranscriptData;

public class GeneTestUtils
{
    public static GeneAnnotation createGeneAnnotation(int svId, boolean isStart, final String geneName, String stableId, int strand,
            final String chromosome, long position, int orientation)
    {
        List<String> synonyms = Lists.newArrayList();
        String karyotypeBand = "";


        GeneAnnotation gene = new GeneAnnotation(svId, isStart, geneName, stableId, strand, synonyms, karyotypeBand);
        gene.setPositionalData(chromosome, position, (byte)orientation);

        return gene;
    }

    // Ensembl data types
    public static EnsemblGeneData createEnsemblGeneData(String geneId, String geneName, String chromosome, int strand, long geneStart, long geneEnd)
    {
        return new EnsemblGeneData(geneId, geneName, chromosome, (byte)strand, geneStart, geneEnd,  "", "");
    }

    public static void addTransExonData(SvGeneTranscriptCollection geneTransCache, final String geneId, List<TranscriptData> transDataList)
    {
        geneTransCache.getTranscriptDataMap().put(geneId, transDataList);
    }

    public static void addGeneData(SvGeneTranscriptCollection geneTransCache, final String chromosome, List<EnsemblGeneData> geneDataList)
    {
        geneTransCache.getChrGeneDataMap().put(chromosome, geneDataList);
    }

    public static TranscriptData createTransExons(final String geneId, int transId, byte strand,
            long[] exonStarts, int[] exonEndPhases, int exonLength)
    {
        return createTransExons(geneId, transId, strand,exonStarts, exonEndPhases, exonLength, false);
    }
    public static TranscriptData createTransExons(final String geneId, int transId, byte strand,
            long[] exonStarts, int[] exonEndPhases, int exonLength, boolean isCanonical)
    {
        if(exonStarts.length == 0 || exonStarts.length != exonEndPhases.length)
            return null;

        String transName = String.format("TRAN%04d", transId);

        int exonCount = exonStarts.length;
        long transStart = exonStarts[0];
        long transEnd = exonStarts[exonCount-1] + exonLength;

        Long codingStart = null;
        Long codingEnd = null;

        int[] exonPhases = new int[exonCount];

        // work out phases and coding start & end
        for(int i = 0; i < exonCount; ++i)
        {
            long exonStart = exonStarts[i];

            int exonEndPhase = exonEndPhases[i];

            if(strand == 1)
                exonPhases[i] = i > 0 ? exonEndPhases[i-1] : -1;
            else
                exonPhases[i] = i < exonCount - 1 ? exonEndPhases[i+1] : -1;

            if(codingStart == null && ((strand == 1 && exonEndPhase != -1) || (strand == -1 && exonPhases[i] != -1)))
            {
                codingStart = new Long(exonStart + exonLength / 2);
            }
            else if(codingStart != null && codingEnd == null
            && ((strand == 1 && exonEndPhase == -1) || (strand == -1 && exonPhases[i] == -1)))
            {
                codingEnd = new Long(exonStart + exonLength / 2);
            }
        }

        TranscriptData transData = new TranscriptData(transId, transName, geneId, isCanonical, strand, transStart, transEnd,
                codingStart, codingEnd, "");

        List<ExonData> exons = Lists.newArrayList();

        for(int i = 0; i < exonCount; ++i)
        {
            long exonStart = exonStarts[i];
            long exonEnd = exonStarts[i] + exonLength;
            int exonRank = strand == 1 ? i + 1 : exonCount - i;

            exons.add(new ExonData(transId, exonStart, exonEnd, exonRank, exonPhases[i], exonEndPhases[i]));
        }

        transData.setExons(exons);

        return transData;

    }

}
