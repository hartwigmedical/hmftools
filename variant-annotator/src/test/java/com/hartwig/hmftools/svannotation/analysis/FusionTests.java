package com.hartwig.hmftools.svannotation.analysis;

import static com.hartwig.hmftools.svannotation.SvGeneTranscriptCollection.extractTranscriptExonData;
import static com.hartwig.hmftools.svannotation.analysis.SvAnnotatorTestUtils.addGeneData;
import static com.hartwig.hmftools.svannotation.analysis.SvAnnotatorTestUtils.addTransExonData;
import static com.hartwig.hmftools.svannotation.analysis.SvAnnotatorTestUtils.createEnsemblGeneData;
import static com.hartwig.hmftools.svannotation.analysis.SvAnnotatorTestUtils.createGeneAnnotation;
import static com.hartwig.hmftools.svannotation.analysis.SvAnnotatorTestUtils.initLogger;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertFalse;
import static org.junit.Assert.assertTrue;

import java.util.List;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.variant.structural.annotation.EnsemblGeneData;
import com.hartwig.hmftools.common.variant.structural.annotation.GeneAnnotation;
import com.hartwig.hmftools.common.variant.structural.annotation.Transcript;
import com.hartwig.hmftools.common.variant.structural.annotation.TranscriptExonData;
import com.hartwig.hmftools.svannotation.SvGeneTranscriptCollection;

import org.junit.Test;

public class FusionTests
{
    @Test
    public void testFusionTypes()
    {
        initLogger();

    }



    @Test
    public void testRnaMatching()
    {
        SvGeneTranscriptCollection geneTransCache = new SvGeneTranscriptCollection();

        String geneName = "GENE1";
        String geneId = "ENSG0001";
        String chromosome = "1";

        List<EnsemblGeneData> geneList = Lists.newArrayList();
        geneList.add(createEnsemblGeneData(geneId, geneName, chromosome, 1, 10000, 20000));

        // one on the negative strand
        String geneName2 = "GENE2";
        String geneId2 = "ENSG0003";
        String chromosome2 = "1";
        geneList.add(createEnsemblGeneData(geneId2, geneName2, chromosome, -1, 10000, 20000));

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


        addTransExonData(geneTransCache, geneId, transExonList);

        String transName2 = "ENST0002";
        int transId2 = 2;

        transExonList = Lists.newArrayList();
        transExonList.add(new TranscriptExonData(geneId2, transName2, transId2, true, (byte)-1, 10100, 19600,
                10500, 10600, 4, -1, -1, null,null, ""));

        transExonList.add(new TranscriptExonData(geneId2, transName2, transId2, true, (byte)-1, 10100, 19600,
                11500, 11600, 3, -1, -1, null, null, ""));

        transExonList.add(new TranscriptExonData(geneId2, transName2, transId2, true, (byte)-1, 10100, 19600,
                12500, 12600, 2, -1, -1, null, null, ""));

        transExonList.add(new TranscriptExonData(geneId2, transName2, transId2, true, (byte)-1, 10100, 19600,
                13500, 13600, 1, -1, -1, null, null, ""));

        addTransExonData(geneTransCache, geneId2, transExonList);


        SvFusionAnalyser fusionAnalyser = new SvFusionAnalyser(null, geneTransCache, "");

        // test positive strand

        long svPos1 = 12700;
        GeneAnnotation geneAnnot1 = createGeneAnnotation(0, true, geneName, geneId, 1, chromosome, svPos1, 1);

        List<TranscriptExonData> transExonDataList = geneTransCache.getTranscriptExons(geneId, transName);
        assertEquals(4, transExonDataList.size());

        Transcript trans = extractTranscriptExonData(transExonDataList, geneAnnot1.position(), geneAnnot1);

        assertTrue(trans != null);

        // test upstream scenarios
        long rnaPosition = 12600;
        boolean isValid = fusionAnalyser.isTranscriptBreakendViableForRnaBoundary(
                trans, true, geneAnnot1.position(), rnaPosition, true);

        assertTrue(isValid);

        // after the next splice site
        svPos1 = 13500;
        geneAnnot1.setPositionalData(chromosome, svPos1, (byte)1);

        isValid = fusionAnalyser.isTranscriptBreakendViableForRnaBoundary(
                trans, true, geneAnnot1.position(), rnaPosition, true);

        assertFalse(isValid);

        // test non-exact RNA boundary
        rnaPosition = 12550;
        isValid = fusionAnalyser.isTranscriptBreakendViableForRnaBoundary(
                trans, true, geneAnnot1.position(), rnaPosition, false);

        assertFalse(isValid);

        rnaPosition = 12700;
        isValid = fusionAnalyser.isTranscriptBreakendViableForRnaBoundary(
                trans, true, geneAnnot1.position(), rnaPosition, false);

        assertFalse(isValid);

        // test downstream

        // exact base at 2nd exon
        svPos1 = 100; // pre promotor
        geneAnnot1.setPositionalData(chromosome, svPos1, (byte)1);

        rnaPosition = 11500;
        isValid = fusionAnalyser.isTranscriptBreakendViableForRnaBoundary(
                trans, false, geneAnnot1.position(), rnaPosition, true);

        assertTrue(isValid);

        // before previous splice acceptor
        svPos1 = 12400;
        geneAnnot1.setPositionalData(chromosome, svPos1, (byte)1);

        rnaPosition = 13500;
        isValid = fusionAnalyser.isTranscriptBreakendViableForRnaBoundary(
                trans, false, geneAnnot1.position(), rnaPosition, true);

        assertFalse(isValid);

        svPos1 = 13000;
        geneAnnot1.setPositionalData(chromosome, svPos1, (byte)1);

        // valid position
        rnaPosition = 13500;
        isValid = fusionAnalyser.isTranscriptBreakendViableForRnaBoundary(
                trans, false, geneAnnot1.position(), rnaPosition, true);

        assertTrue(isValid);


        // now test the negative strand

        long svPos2 = 12700;
        GeneAnnotation geneAnnot2 = createGeneAnnotation(1, true, geneName2, geneId2, -1, chromosome, svPos2, -1);

        transExonDataList = geneTransCache.getTranscriptExons(geneId2, transName2);
        assertEquals(4, transExonDataList.size());

        Transcript trans2 = extractTranscriptExonData(transExonDataList, geneAnnot2.position(), geneAnnot2);

        assertTrue(trans2 != null);

        // upstream

        rnaPosition = 11500; // 3rd exon end

        svPos2 = 11600;
        geneAnnot2.setPositionalData(chromosome, svPos2, (byte)1);

        isValid = fusionAnalyser.isTranscriptBreakendViableForRnaBoundary(
                trans2, true, geneAnnot2.position(), rnaPosition, true);

        assertTrue(isValid);

        // test downstream

        rnaPosition = 11600; // 3rd exon start

        svPos2 = 11700;
        geneAnnot2.setPositionalData(chromosome, svPos2, (byte)1);

        isValid = fusionAnalyser.isTranscriptBreakendViableForRnaBoundary(
                trans2, false, geneAnnot2.position(), rnaPosition, true);

        assertTrue(isValid);

        // before prev splice acceptor is invali
        svPos2 = 12700;
        geneAnnot2.setPositionalData(chromosome, svPos2, (byte)1);

        isValid = fusionAnalyser.isTranscriptBreakendViableForRnaBoundary(
                trans2, false, geneAnnot2.position(), rnaPosition, true);

        assertFalse(isValid);

        // invalid too far upstream of promotor
        rnaPosition = 12600; // 3rd exon start

        svPos2 = 130000;
        geneAnnot2.setPositionalData(chromosome, svPos2, (byte)1);

        isValid = fusionAnalyser.isTranscriptBreakendViableForRnaBoundary(
                trans2, false, geneAnnot2.position(), rnaPosition, true);

        assertFalse(isValid);

    }
}
