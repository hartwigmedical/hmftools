package com.hartwig.hmftools.svanalysis.analyser.com.hartwig.hmftools.svanalysis.gene;

import static com.hartwig.hmftools.svanalysis.analyser.com.hartwig.hmftools.svanalysis.gene.GeneTestUtils.addGeneData;
import static com.hartwig.hmftools.svanalysis.analyser.com.hartwig.hmftools.svanalysis.gene.GeneTestUtils.addTransExonData;
import static com.hartwig.hmftools.svanalysis.analyser.com.hartwig.hmftools.svanalysis.gene.GeneTestUtils.createEnsemblGeneData;
import static com.hartwig.hmftools.svanalysis.analyser.com.hartwig.hmftools.svanalysis.gene.GeneTestUtils.createGeneAnnotation;
import static com.hartwig.hmftools.svanalysis.gene.SvGeneTranscriptCollection.EXON_RANK_MAX;
import static com.hartwig.hmftools.svanalysis.gene.SvGeneTranscriptCollection.EXON_RANK_MIN;
import static com.hartwig.hmftools.svanalysis.gene.SvGeneTranscriptCollection.extractTranscriptExonData;
import static com.hartwig.hmftools.svanalysis.gene.SvGeneTranscriptCollection.setAlternativeTranscriptPhasings;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertTrue;

import java.util.List;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.variant.structural.annotation.EnsemblGeneData;
import com.hartwig.hmftools.common.variant.structural.annotation.GeneAnnotation;
import com.hartwig.hmftools.common.variant.structural.annotation.Transcript;
import com.hartwig.hmftools.common.variant.structural.annotation.TranscriptExonData;
import com.hartwig.hmftools.svanalysis.gene.SvGeneTranscriptCollection;

import org.junit.Test;

public class GeneCollectionTest
{
    @Test
    public void testExonDataExtraction()
    {
        SvGeneTranscriptCollection geneTransCache = new SvGeneTranscriptCollection();

        String geneName = "GENE1";
        String geneId = "ENSG0001";
        String chromosome = "1";

        List<EnsemblGeneData> geneList = Lists.newArrayList();
        geneList.add(createEnsemblGeneData(geneId, geneName, chromosome, 1, 10000, 20000));
        addGeneData(geneTransCache, chromosome, geneList);

        List<TranscriptExonData> transExonList = Lists.newArrayList();

        String transName = "ENST0001";
        int transId = 1;

        long codingStart = 11550;
        long codingEnd = 13550;
        transExonList.add(new TranscriptExonData(geneId, transName, transId, true, (byte)1, 10100, 19600,
                10500, 10600, 1, -1, -1, codingStart, codingEnd, ""));

        transExonList.add(new TranscriptExonData(geneId, transName, transId, true, (byte)1, 10100, 19600,
                11500, 11600, 2, -1, 1, codingStart, codingEnd, ""));

        transExonList.add(new TranscriptExonData(geneId, transName, transId, true, (byte)1, 10100, 19600,
                12500, 12600, 3, 1, 2, codingStart, codingEnd, ""));

        transExonList.add(new TranscriptExonData(geneId, transName, transId, true, (byte)1, 10100, 19600,
                13500, 13600, 4, 2, -1, codingStart, codingEnd, ""));

        String transName2 = "ENST0002";
        int transId2 = 2;

        transExonList.add(new TranscriptExonData(geneId, transName2, transId2, false, (byte)1, 10100, 19600,
                12500, 12600, 1, -1, -1, null,null, ""));

        transExonList.add(new TranscriptExonData(geneId, transName2, transId2, false, (byte)1, 10100, 19600,
                13500, 13600, 2, -1, -1, null, null, ""));

        transExonList.add(new TranscriptExonData(geneId, transName2, transId2, false, (byte)1, 10100, 19600,
                14500, 14600, 3, -1, -1, null, null, ""));

        addTransExonData(geneTransCache, geneId, transExonList);

        // test exon retrieval
        List<TranscriptExonData> exonDataList = geneTransCache.getTranscriptExons(geneId, "");
        assertEquals(4, exonDataList.size());
        assertEquals(transId, exonDataList.get(0).TransId);

        exonDataList = geneTransCache.getTranscriptExons(geneId, transName2);
        assertEquals(3, exonDataList.size());
        assertEquals(transId2, exonDataList.get(0).TransId);

        // test exon ranks given a position

        int[] transUpExonData = geneTransCache.getExonRankings(geneId, 11400);

        assertEquals(1, transUpExonData[EXON_RANK_MIN]);
        assertEquals(2, transUpExonData[EXON_RANK_MAX]);

        // before the first
        transUpExonData = geneTransCache.getExonRankings(geneId, 9000);

        assertEquals(0, transUpExonData[EXON_RANK_MIN]);
        assertEquals(1, transUpExonData[EXON_RANK_MAX]);

        // after the last
        transUpExonData = geneTransCache.getExonRankings(geneId, 16000);

        assertEquals(4, transUpExonData[EXON_RANK_MIN]);
        assertEquals(-1, transUpExonData[EXON_RANK_MAX]);


        // on an exon boundary
        transUpExonData = geneTransCache.getExonRankings(geneId, 12500);

        assertEquals(3, transUpExonData[EXON_RANK_MIN]);
        assertEquals(3, transUpExonData[EXON_RANK_MAX]);

    }

    @Test
    public void testTranscriptBreakends()
    {
        SvGeneTranscriptCollection geneTransCache = new SvGeneTranscriptCollection();

        // first a gene on the forward strand
        String geneName = "GENE1";
        String geneId = "ENSG0001";
        String chromosome = "1";

        List<EnsemblGeneData> geneList = Lists.newArrayList();
        geneList.add(createEnsemblGeneData(geneId, geneName, chromosome, 1, 100, 1000));
        addGeneData(geneTransCache, chromosome, geneList);

        GeneAnnotation genePosStrand = createGeneAnnotation(0, true, geneName, geneId, 1, chromosome, 0, 1);

        List<TranscriptExonData> transExonList = Lists.newArrayList();

        String transName = "ENST0001";
        int transId = 1;

        long codingStart = 350;
        long codingEnd = 950;
        long transStart = 100;
        long transEnd = 1000;
        transExonList.add(new TranscriptExonData(geneId, transName, transId, true, (byte)1, transStart, transEnd,
                100, 200, 1, -1, -1, codingStart, codingEnd, ""));

        transExonList.add(new TranscriptExonData(geneId, transName, transId, true, (byte)1, transStart, transEnd,
                300, 400, 2, -1, 1, codingStart, codingEnd, ""));

        transExonList.add(new TranscriptExonData(geneId, transName, transId, true, (byte)1, transStart, transEnd,
                500, 600, 3, 1, 2, codingStart, codingEnd, ""));

        transExonList.add(new TranscriptExonData(geneId, transName, transId, true, (byte)1, transStart, transEnd,
                700, 800, 4, 2, 0, codingStart, codingEnd, ""));

        transExonList.add(new TranscriptExonData(geneId, transName, transId, true, (byte)1, transStart, transEnd,
                900, 1000, 5, 0, -1, codingStart, codingEnd, ""));

        addTransExonData(geneTransCache, geneId, transExonList);

        long position = 250;
        byte posOrientation = 1;
        byte negOrientation = -1;
        Transcript trans = extractTranscriptExonData(transExonList, position, genePosStrand);

        assertEquals(5, trans.ExonMax);
        assertEquals(1, trans.ExonUpstream);
        assertEquals(2, trans.ExonDownstream);
        assertEquals(-1, trans.ExonUpstreamPhase);
        assertEquals(-1, trans.ExonDownstreamPhase);

        // test caching of upstream phasings for exon-skipping fusion logic
        setAlternativeTranscriptPhasings(trans, transExonList, position, posOrientation);
        assertEquals(0, trans.getAlternativePhasing().size());

        // and test as a downstream gene
        setAlternativeTranscriptPhasings(trans, transExonList, position, negOrientation);
        assertEquals(3, trans.getAlternativePhasing().size());
        Integer exonsSkipped = trans.getAlternativePhasing().get(1);
        assertTrue(exonsSkipped != null && exonsSkipped == 1);
        exonsSkipped = trans.getAlternativePhasing().get(2);
        assertTrue(exonsSkipped != null && exonsSkipped == 2);
        exonsSkipped = trans.getAlternativePhasing().get(0);
        assertTrue(exonsSkipped != null && exonsSkipped == 3);

        position = 450;
        trans = extractTranscriptExonData(transExonList, position, genePosStrand);

        assertEquals(2, trans.ExonUpstream);
        assertEquals(3, trans.ExonDownstream);
        assertEquals(1, trans.ExonUpstreamPhase);
        assertEquals(1, trans.ExonDownstreamPhase);

        setAlternativeTranscriptPhasings(trans, transExonList, position, posOrientation);
        assertEquals(1, trans.getAlternativePhasing().size());
        exonsSkipped = trans.getAlternativePhasing().get(-1);
        assertTrue(exonsSkipped != null && exonsSkipped == 1);

        setAlternativeTranscriptPhasings(trans, transExonList, position, negOrientation);
        assertEquals(2, trans.getAlternativePhasing().size());
        exonsSkipped = trans.getAlternativePhasing().get(2);
        assertTrue(exonsSkipped != null && exonsSkipped == 1);
        exonsSkipped = trans.getAlternativePhasing().get(0);
        assertTrue(exonsSkipped != null && exonsSkipped == 2);

        position = 650;
        trans = extractTranscriptExonData(transExonList, position, genePosStrand);

        assertEquals(3, trans.ExonUpstream);
        assertEquals(4, trans.ExonDownstream);
        assertEquals(2, trans.ExonUpstreamPhase);
        assertEquals(2, trans.ExonDownstreamPhase);

        setAlternativeTranscriptPhasings(trans, transExonList, position, posOrientation);
        assertEquals(2, trans.getAlternativePhasing().size());
        exonsSkipped = trans.getAlternativePhasing().get(-1);
        assertTrue(exonsSkipped != null && exonsSkipped == 2);
        exonsSkipped = trans.getAlternativePhasing().get(1);
        assertTrue(exonsSkipped != null && exonsSkipped == 1);

        setAlternativeTranscriptPhasings(trans, transExonList, position, negOrientation);
        assertEquals(1, trans.getAlternativePhasing().size());
        exonsSkipped = trans.getAlternativePhasing().get(0);
        assertTrue(exonsSkipped != null && exonsSkipped == 1);


        // then a gene on the reverse strand
        geneName = "GENE2";
        geneId = "ENSG0002";
        chromosome = "1";

        geneList.add(createEnsemblGeneData(geneId, geneName, chromosome, 1, 100, 1000));
        addGeneData(geneTransCache, chromosome, geneList);

        GeneAnnotation geneNegStrand = createGeneAnnotation(0, true, geneName, geneId, -1, chromosome, 0, 1);

        transExonList = Lists.newArrayList();

        transName = "ENST0002";
        transId = 2;

        codingStart = 350;
        codingEnd = 950;

        transExonList.add(new TranscriptExonData(geneId, transName, transId, true, (byte)1, transStart, transEnd,
                100, 200, 5, -1, -1, codingStart, codingEnd, ""));

        transExonList.add(new TranscriptExonData(geneId, transName, transId, true, (byte)1, transStart, transEnd,
                300, 400, 4, 1, -1, codingStart, codingEnd, ""));

        transExonList.add(new TranscriptExonData(geneId, transName, transId, true, (byte)1, transStart, transEnd,
                500, 600, 3, 0, 1, codingStart, codingEnd, ""));

        transExonList.add(new TranscriptExonData(geneId, transName, transId, true, (byte)1, transStart, transEnd,
                700, 800, 2, 2, 0, codingStart, codingEnd, ""));

        transExonList.add(new TranscriptExonData(geneId, transName, transId, true, (byte)1, transStart, transEnd,
                900, 1000, 1, -1, 2, codingStart, codingEnd, ""));

        addTransExonData(geneTransCache, geneId, transExonList);

        position = 850;
        trans = extractTranscriptExonData(transExonList, position, geneNegStrand);

        assertEquals(5, trans.ExonMax);
        assertEquals(1, trans.ExonUpstream);
        assertEquals(2, trans.ExonDownstream);
        assertEquals(2, trans.ExonUpstreamPhase);
        assertEquals(2, trans.ExonDownstreamPhase);

        setAlternativeTranscriptPhasings(trans, transExonList, position, negOrientation);
        assertEquals(0, trans.getAlternativePhasing().size());

        setAlternativeTranscriptPhasings(trans, transExonList, position, posOrientation);
        assertEquals(3, trans.getAlternativePhasing().size());
        exonsSkipped = trans.getAlternativePhasing().get(0);
        assertTrue(exonsSkipped != null && exonsSkipped == 1);
        exonsSkipped = trans.getAlternativePhasing().get(1);
        assertTrue(exonsSkipped != null && exonsSkipped == 2);
        exonsSkipped = trans.getAlternativePhasing().get(-1);
        assertTrue(exonsSkipped != null && exonsSkipped == 3);

        position = 250;
        trans = extractTranscriptExonData(transExonList, position, geneNegStrand);

        assertEquals(4, trans.ExonUpstream);
        assertEquals(5, trans.ExonDownstream);
        assertEquals(-1, trans.ExonUpstreamPhase);
        assertEquals(-1, trans.ExonDownstreamPhase);

        setAlternativeTranscriptPhasings(trans, transExonList, position, negOrientation);
        assertEquals(3, trans.getAlternativePhasing().size());
        exonsSkipped = trans.getAlternativePhasing().get(1);
        assertTrue(exonsSkipped != null && exonsSkipped == 1);
        exonsSkipped = trans.getAlternativePhasing().get(0);
        assertTrue(exonsSkipped != null && exonsSkipped == 2);
        exonsSkipped = trans.getAlternativePhasing().get(2);
        assertTrue(exonsSkipped != null && exonsSkipped == 3);

        setAlternativeTranscriptPhasings(trans, transExonList, position, posOrientation);
        assertEquals(0, trans.getAlternativePhasing().size());
    }

}
