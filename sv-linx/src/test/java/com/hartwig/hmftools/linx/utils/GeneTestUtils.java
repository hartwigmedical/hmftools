package com.hartwig.hmftools.linx.utils;

import static com.hartwig.hmftools.common.drivercatalog.DriverType.DEL;
import static com.hartwig.hmftools.common.drivercatalog.LikelihoodMethod.AMP;
import static com.hartwig.hmftools.common.purple.copynumber.CopyNumberMethod.BAF_WEIGHTED;

import java.util.List;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.drivercatalog.DriverCatalog;
import com.hartwig.hmftools.common.drivercatalog.DriverCategory;
import com.hartwig.hmftools.common.drivercatalog.DriverType;
import com.hartwig.hmftools.common.drivercatalog.ImmutableDriverCatalog;
import com.hartwig.hmftools.common.drivercatalog.LikelihoodMethod;
import com.hartwig.hmftools.common.purple.gene.GeneCopyNumber;
import com.hartwig.hmftools.common.purple.gene.ImmutableGeneCopyNumber;
import com.hartwig.hmftools.common.purple.segment.SegmentSupport;
import com.hartwig.hmftools.common.variant.structural.annotation.EnsemblGeneData;
import com.hartwig.hmftools.common.variant.structural.annotation.ExonData;
import com.hartwig.hmftools.common.variant.structural.annotation.GeneAnnotation;
import com.hartwig.hmftools.common.variant.structural.annotation.TranscriptData;
import com.hartwig.hmftools.linx.gene.SvGeneTranscriptCollection;

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

    public static String generateTransName(int transId) { return String.format("TRAN%04d", transId); }

    public static TranscriptData createTransExons(final String geneId, int transId, byte strand,
            long[] exonStarts, int[] exonEndPhases, int exonLength, boolean isCanonical)
    {
        if(exonStarts.length == 0 || exonStarts.length != exonEndPhases.length)
            return null;

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

        TranscriptData transData = new TranscriptData(transId, generateTransName(transId), geneId, isCanonical, strand, transStart, transEnd,
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

        public static GeneCopyNumber createGeneCopyNumber(final String gene, final String chromosome,
                double minCopyNumber, long posStart, long posEnd)
    {
        return ImmutableGeneCopyNumber.builder()
                .chromosome(chromosome)
                .chromosomeBand("")
                .start(posStart)
                .end(posEnd)
                .gene(gene)
                .maxCopyNumber(minCopyNumber)
                .minCopyNumber(minCopyNumber)
                .somaticRegions(0)
                .germlineHet2HomRegions(0)
                .germlineHomRegions(0)
                .minRegions(0)
                .minRegionStart(0)
                .minRegionEnd(0)
                .minRegionStartSupport(SegmentSupport.BND)
                .minRegionEndSupport(SegmentSupport.BND)
                .minRegionMethod(BAF_WEIGHTED)
                .minMinorAllelePloidy(0)
                .transcriptID("")
                .transcriptVersion(1)
                .build();
    }

    public static DriverCatalog createDriver(final String gene, final String chromosome, DriverType type, DriverCategory category,
            boolean biallelic, double minCopyNumber)
    {
        LikelihoodMethod method;
        if(type == DriverType.AMP)
            method = LikelihoodMethod.AMP;
        else if(type == DriverType.DEL)
            method = LikelihoodMethod.DEL;
        else
            method = LikelihoodMethod.BIALLELIC;

        return ImmutableDriverCatalog.builder()
                .biallelic(biallelic)
                .category(category)
                .gene(gene)
                .chromosome(chromosome)
                .chromosomeBand("")
                .driver(type)
                .driverLikelihood(1.0)
                .dndsLikelihood(1.0)
                .likelihoodMethod(method)
                .minCopyNumber(minCopyNumber)
                .maxCopyNumber(minCopyNumber)
                .missense(0)
                .nonsense(0)
                .splice(0)
                .inframe(0)
                .frameshift(0)
                .build();
    }

}
