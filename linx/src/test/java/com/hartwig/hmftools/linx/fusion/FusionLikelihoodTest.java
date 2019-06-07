package com.hartwig.hmftools.linx.fusion;

import static com.hartwig.hmftools.linx.analyser.gene.GeneTestUtils.createEnsemblGeneData;
import static com.hartwig.hmftools.linx.fusion.FusionLikelihood.getGenePhasingCounts;
import static com.hartwig.hmftools.linx.fusion.GeneRangeData.GENE_PHASING_REGION_5P_UTR;

import static org.junit.Assert.assertEquals;

import java.util.List;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.variant.structural.annotation.EnsemblGeneData;
import com.hartwig.hmftools.common.variant.structural.annotation.TranscriptExonData;

import org.junit.Test;

public class FusionLikelihoodTest
{
    @Test
    public void testGeneRegionCounts()
    {
        String geneId  = "G001";
        byte strand = 1;
        EnsemblGeneData geneData = createEnsemblGeneData(geneId, "GEN1", "1", strand, 10, 100);

        List<TranscriptExonData> transExonDataList = Lists.newArrayList();

        // 3 coding exons, with coding region starting and ending half way through them

        int transStart = 10;
        int transEnd = 100;
        Long codingStart = new Long(35);
        Long codingEnd = new Long(75);

        transExonDataList.add(new TranscriptExonData(geneId, "T001", 1, true, strand, transStart, transEnd,
                10, 20, 1, -1, -1, codingStart, codingEnd, ""));

        transExonDataList.add(new TranscriptExonData(geneId, "T001", 1, true, strand, transStart, transEnd,
                30, 40, 1, -1, 1, codingStart, codingEnd, ""));

        transExonDataList.add(new TranscriptExonData(geneId, "T001", 1, true, strand, transStart, transEnd,
                50, 60, 1, 1, 2, codingStart, codingEnd, ""));

        transExonDataList.add(new TranscriptExonData(geneId, "T001", 1, true, strand, transStart, transEnd,
                70, 80, 1, 2, -1, codingStart, codingEnd, ""));

        transExonDataList.add(new TranscriptExonData(geneId, "T001", 1, true, strand, transStart, transEnd,
                90, 100, 1, -1, -1, codingStart, codingEnd, ""));

        int[] regionTotals = getGenePhasingCounts(geneData, transExonDataList);

        assertEquals(regionTotals[GENE_PHASING_REGION_5P_UTR], 26);
        //assertEquals(regionTotals[GENE_PHASING_REGION_CODING_0], 7);
        //assertEquals(regionTotals[GENE_PHASING_REGION_CODING_1], 17);
        //assertEquals(regionTotals[GENE_PHASING_REGION_CODING_2], 16);


        // now test with the reverse strand
        geneId  = "G002";
        strand = -1;
        EnsemblGeneData geneData2 = createEnsemblGeneData(geneId, "GEN2", "1", strand, 10, 100);

        transExonDataList.clear();

        // 3 coding exons, with coding region starting and ending half way through them

        codingStart = new Long(35);
        codingEnd = new Long(75);

        String transId = "T001";

        transExonDataList.add(new TranscriptExonData(geneId, transId, 1, true, strand, transStart, transEnd,
                10, 20, 1, -1, -1, codingStart, codingEnd, ""));

        transExonDataList.add(new TranscriptExonData(geneId, transId, 1, true, strand, transStart, transEnd,
                30, 40, 1, -1, 1, codingStart, codingEnd, ""));

        transExonDataList.add(new TranscriptExonData(geneId, transId, 1, true, strand, transStart, transEnd,
                50, 60, 1, 1, 2, codingStart, codingEnd, ""));

        transExonDataList.add(new TranscriptExonData(geneId, transId, 1, true, strand, transStart, transEnd,
                70, 80, 1, 2, -1, codingStart, codingEnd, ""));

        transExonDataList.add(new TranscriptExonData(geneId, transId, 1, true, strand, transStart, transEnd,
                90, 100, 1, -1, -1, codingStart, codingEnd, ""));

        regionTotals = getGenePhasingCounts(geneData2, transExonDataList);

        // assertEquals(regionTotals[GENE_PHASING_REGION_5P_UTR], 26);
        //assertEquals(regionTotals[GENE_PHASING_REGION_CODING_0], 7);
        //assertEquals(regionTotals[GENE_PHASING_REGION_CODING_1], 17);
        //assertEquals(regionTotals[GENE_PHASING_REGION_CODING_2], 16);

    }

}
