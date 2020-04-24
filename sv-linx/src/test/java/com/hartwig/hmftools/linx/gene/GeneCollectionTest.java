package com.hartwig.hmftools.linx.gene;

import static com.hartwig.hmftools.common.ensemblcache.EnsemblDataCache.EXON_PHASE_MAX;
import static com.hartwig.hmftools.common.ensemblcache.EnsemblDataCache.EXON_PHASE_MIN;
import static com.hartwig.hmftools.common.ensemblcache.EnsemblDataCache.EXON_RANK_MAX;
import static com.hartwig.hmftools.common.ensemblcache.EnsemblDataCache.EXON_RANK_MIN;
import static com.hartwig.hmftools.common.ensemblcache.EnsemblDataCache.extractTranscriptExonData;
import static com.hartwig.hmftools.common.ensemblcache.EnsemblDataCache.getProteinDomainPositions;
import static com.hartwig.hmftools.common.ensemblcache.EnsemblDataCache.setAlternativeTranscriptPhasings;
import static com.hartwig.hmftools.common.ensemblcache.TranscriptProteinData.BIOTYPE_PROTEIN_CODING;
import static com.hartwig.hmftools.common.variant.structural.annotation.Transcript.POST_CODING_PHASE;
import static com.hartwig.hmftools.common.variant.structural.annotation.Transcript.TRANS_CODING_TYPE_3P_UTR;
import static com.hartwig.hmftools.common.variant.structural.annotation.Transcript.TRANS_CODING_TYPE_CODING;
import static com.hartwig.hmftools.common.ensemblcache.GeneTestUtils.createGeneDataCache;
import static com.hartwig.hmftools.common.ensemblcache.GeneTestUtils.createTransExons;
import static com.hartwig.hmftools.common.utils.sv.StartEndIterator.SE_END;
import static com.hartwig.hmftools.common.utils.sv.StartEndIterator.SE_START;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertTrue;

import java.util.List;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.ensemblcache.EnsemblDataCache;
import com.hartwig.hmftools.common.ensemblcache.EnsemblGeneData;
import com.hartwig.hmftools.common.ensemblcache.GeneTestUtils;
import com.hartwig.hmftools.common.variant.structural.annotation.GeneAnnotation;
import com.hartwig.hmftools.common.variant.structural.annotation.Transcript;
import com.hartwig.hmftools.common.ensemblcache.TranscriptData;
import com.hartwig.hmftools.common.ensemblcache.TranscriptProteinData;

import org.junit.Test;

public class GeneCollectionTest
{
    @Test
    public void testExonDataExtraction()
    {
        EnsemblDataCache geneTransCache = createGeneDataCache();

        String geneName = "GENE1";
        String geneId = "ENSG0001";
        String chromosome = "1";

        List<EnsemblGeneData> geneList = Lists.newArrayList();
        geneList.add(GeneTestUtils.createEnsemblGeneData(geneId, geneName, chromosome, 1, 10000, 20000));
        GeneTestUtils.addGeneData(geneTransCache, chromosome, geneList);

        List<TranscriptData> transDataList = Lists.newArrayList();

        int transId = 1;
        byte strand = 1;

        int[] exonStarts = new int[]{10500, 11500, 12500, 13500};
        int[] exonPhases = new int[]{-1, 1, 2, -1};

        TranscriptData transData = createTransExons(geneId, transId, strand, exonStarts, exonPhases, 100, true);
        transDataList.add(transData);

        int transId2 = 2;

        exonStarts = new int[]{12500, 13500, 14500};
        exonPhases = new int[]{-1, -1, -1};

        transData = createTransExons(geneId, transId2, strand, exonStarts, exonPhases, 100);
        String transName2 = transData.TransName;
        transDataList.add(transData);

        GeneTestUtils.addTransExonData(geneTransCache, geneId, transDataList);

        // test exon retrieval
        transData = geneTransCache.getTranscriptData(geneId, "");
        assertEquals(transId, transData.TransId);
        assertEquals(4, transData.exons().size());

        transData = geneTransCache.getTranscriptData(geneId, transName2);
        assertEquals(transId2, transData.TransId);
        assertEquals(3, transData.exons().size());

        // test exon ranks given a position

        int[] transUpExonData = geneTransCache.getExonRankings(geneId, 11400);

        assertEquals(1, transUpExonData[EXON_RANK_MIN]);
        assertEquals(2, transUpExonData[EXON_RANK_MAX]);
        assertEquals(-1, transUpExonData[EXON_PHASE_MIN]);
        assertEquals(-1, transUpExonData[EXON_PHASE_MAX]);

        // before the first
        transUpExonData = geneTransCache.getExonRankings(geneId, 9000);

        assertEquals(0, transUpExonData[EXON_RANK_MIN]);
        assertEquals(1, transUpExonData[EXON_RANK_MAX]);
        assertEquals(-1, transUpExonData[EXON_PHASE_MIN]);
        assertEquals(-1, transUpExonData[EXON_PHASE_MAX]);

        // after the last
        transUpExonData = geneTransCache.getExonRankings(geneId, 16000);

        assertEquals(4, transUpExonData[EXON_RANK_MIN]);
        assertEquals(-1, transUpExonData[EXON_RANK_MAX]);
        assertEquals(-1, transUpExonData[EXON_PHASE_MIN]);
        assertEquals(-1, transUpExonData[EXON_PHASE_MAX]);

        // on an exon boundary
        transUpExonData = geneTransCache.getExonRankings(geneId, 12500);

        assertEquals(3, transUpExonData[EXON_RANK_MIN]);
        assertEquals(3, transUpExonData[EXON_RANK_MAX]);
        assertEquals(1, transUpExonData[EXON_PHASE_MIN]);
        assertEquals(1, transUpExonData[EXON_PHASE_MAX]);
    }

    @Test
    public void testTranscriptBreakends()
    {
        EnsemblDataCache geneTransCache = createGeneDataCache();

        // first a gene on the forward strand
        String geneName = "GENE1";
        String geneId = "ENSG0001";
        String chromosome = "1";

        List<EnsemblGeneData> geneList = Lists.newArrayList();
        geneList.add(GeneTestUtils.createEnsemblGeneData(geneId, geneName, chromosome, 1, 100, 1000));
        GeneTestUtils.addGeneData(geneTransCache, chromosome, geneList);

        GeneAnnotation genePosStrand = GeneTestUtils.createGeneAnnotation(0, true, geneName, geneId, 1, chromosome, 0, 1);

        List<TranscriptData> transDataList = Lists.newArrayList();

        int transId = 1;
        byte strand = 1;

        int[] exonStarts = new int[]{100, 300, 500, 700, 900};
        int[] exonPhases = new int[]{-1, 1, 2, 0, -1};

        TranscriptData transData = createTransExons(geneId, transId++, strand, exonStarts, exonPhases, 100);
        transDataList.add(transData);

        GeneTestUtils.addTransExonData(geneTransCache, geneId, transDataList);

        int position = 250;
        byte posOrientation = 1;
        byte negOrientation = -1;
        Transcript trans = extractTranscriptExonData(transData, position, genePosStrand);

        assertEquals(5, trans.ExonMax);
        assertEquals(1, trans.ExonUpstream);
        assertEquals(2, trans.ExonDownstream);
        assertEquals(-1, trans.ExonUpstreamPhase);
        assertEquals(-1, trans.ExonDownstreamPhase);

        // test caching of upstream phasings for exon-skipping fusion logic
        setAlternativeTranscriptPhasings(trans, transData.exons(), position, posOrientation);
        assertEquals(0, trans.getAlternativePhasing().size());

        // and test as a downstream gene
        setAlternativeTranscriptPhasings(trans, transData.exons(), position, negOrientation);
        assertEquals(3, trans.getAlternativePhasing().size());
        Integer exonsSkipped = trans.getAlternativePhasing().get(1);
        assertTrue(exonsSkipped != null && exonsSkipped == 1);
        exonsSkipped = trans.getAlternativePhasing().get(2);
        assertTrue(exonsSkipped != null && exonsSkipped == 2);
        exonsSkipped = trans.getAlternativePhasing().get(0);
        assertTrue(exonsSkipped != null && exonsSkipped == 3);

        position = 450;
        trans = extractTranscriptExonData(transData, position, genePosStrand);

        assertEquals(2, trans.ExonUpstream);
        assertEquals(3, trans.ExonDownstream);
        assertEquals(1, trans.ExonUpstreamPhase);
        assertEquals(1, trans.ExonDownstreamPhase);

        setAlternativeTranscriptPhasings(trans, transData.exons(), position, posOrientation);
        assertEquals(1, trans.getAlternativePhasing().size());
        exonsSkipped = trans.getAlternativePhasing().get(-1);
        assertTrue(exonsSkipped != null && exonsSkipped == 1);

        setAlternativeTranscriptPhasings(trans, transData.exons(), position, negOrientation);
        assertEquals(2, trans.getAlternativePhasing().size());
        exonsSkipped = trans.getAlternativePhasing().get(2);
        assertTrue(exonsSkipped != null && exonsSkipped == 1);
        exonsSkipped = trans.getAlternativePhasing().get(0);
        assertTrue(exonsSkipped != null && exonsSkipped == 2);

        position = 650;
        trans = extractTranscriptExonData(transData, position, genePosStrand);

        assertEquals(3, trans.ExonUpstream);
        assertEquals(4, trans.ExonDownstream);
        assertEquals(2, trans.ExonUpstreamPhase);
        assertEquals(2, trans.ExonDownstreamPhase);

        setAlternativeTranscriptPhasings(trans, transData.exons(), position, posOrientation);
        assertEquals(2, trans.getAlternativePhasing().size());
        exonsSkipped = trans.getAlternativePhasing().get(-1);
        assertTrue(exonsSkipped != null && exonsSkipped == 2);
        exonsSkipped = trans.getAlternativePhasing().get(1);
        assertTrue(exonsSkipped != null && exonsSkipped == 1);

        setAlternativeTranscriptPhasings(trans, transData.exons(), position, negOrientation);
        assertEquals(1, trans.getAlternativePhasing().size());
        exonsSkipped = trans.getAlternativePhasing().get(0);
        assertTrue(exonsSkipped != null && exonsSkipped == 1);


        // then a gene on the reverse strand
        geneName = "GENE2";
        geneId = "ENSG0002";
        chromosome = "1";

        geneList.add(GeneTestUtils.createEnsemblGeneData(geneId, geneName, chromosome, 1, 100, 1000));
        GeneTestUtils.addGeneData(geneTransCache, chromosome, geneList);

        GeneAnnotation geneNegStrand = GeneTestUtils.createGeneAnnotation(0, true, geneName, geneId, -1, chromosome, 0, 1);

        transDataList = Lists.newArrayList();

        transId = 2;
        strand = -1;

        exonStarts = new int[]{100, 300, 500, 700, 900};
        exonPhases = new int[]{-1, -1, 1, 0, 2};

        transData = createTransExons(geneId, transId++, strand, exonStarts, exonPhases, 100);
        transDataList.add(transData);

        GeneTestUtils.addTransExonData(geneTransCache, geneId, transDataList);

        position = 850;
        trans = extractTranscriptExonData(transData, position, geneNegStrand);

        assertEquals(5, trans.ExonMax);
        assertEquals(1, trans.ExonUpstream);
        assertEquals(2, trans.ExonDownstream);
        assertEquals(2, trans.ExonUpstreamPhase);
        assertEquals(2, trans.ExonDownstreamPhase);

        setAlternativeTranscriptPhasings(trans, transData.exons(), position, negOrientation);
        assertEquals(0, trans.getAlternativePhasing().size());

        setAlternativeTranscriptPhasings(trans, transData.exons(), position, posOrientation);
        assertEquals(3, trans.getAlternativePhasing().size());
        exonsSkipped = trans.getAlternativePhasing().get(0);
        assertTrue(exonsSkipped != null && exonsSkipped == 1);
        exonsSkipped = trans.getAlternativePhasing().get(1);
        assertTrue(exonsSkipped != null && exonsSkipped == 2);
        exonsSkipped = trans.getAlternativePhasing().get(-1);
        assertTrue(exonsSkipped != null && exonsSkipped == 3);

        position = 250;
        trans = extractTranscriptExonData(transData, position, geneNegStrand);

        assertEquals(4, trans.ExonUpstream);
        assertEquals(5, trans.ExonDownstream);
        assertEquals(POST_CODING_PHASE, trans.ExonUpstreamPhase);
        assertEquals(POST_CODING_PHASE, trans.ExonDownstreamPhase);

        setAlternativeTranscriptPhasings(trans, transData.exons(), position, negOrientation);
        assertEquals(3, trans.getAlternativePhasing().size());
        exonsSkipped = trans.getAlternativePhasing().get(1);
        assertTrue(exonsSkipped != null && exonsSkipped == 1);
        exonsSkipped = trans.getAlternativePhasing().get(0);
        assertTrue(exonsSkipped != null && exonsSkipped == 2);
        exonsSkipped = trans.getAlternativePhasing().get(2);
        assertTrue(exonsSkipped != null && exonsSkipped == 3);

        setAlternativeTranscriptPhasings(trans, transData.exons(), position, posOrientation);
        assertEquals(0, trans.getAlternativePhasing().size());
    }

    @Test
    public void testBreakendTranscriptCoding()
    {
        EnsemblDataCache geneTransCache = createGeneDataCache();

        // first a gene on the forward strand
        String geneName = "GENE1";
        String geneId = "ENSG0001";
        String chromosome = "1";

        GeneAnnotation genePosStrand = GeneTestUtils.createGeneAnnotation(0, true, geneName, geneId, 1, chromosome, 0, 1);

        int transId = 1;
        byte strand = 1;

        int[] exonStarts = new int[]{100, 200, 300, 400, 500};

        // coding taking up exactly the first exon
        Integer codingStart = new Integer(100);
        Integer codingEnd = new Integer(110);

        TranscriptData transData = createTransExons(
                geneId, transId++, strand, exonStarts, 10, codingStart, codingEnd, true, BIOTYPE_PROTEIN_CODING);

        int position = 150;
        Transcript trans = extractTranscriptExonData(transData, position, genePosStrand);

        assertEquals(5, trans.ExonMax);
        assertEquals(1, trans.ExonUpstream);
        assertEquals(2, trans.ExonDownstream);
        assertEquals(TRANS_CODING_TYPE_3P_UTR, trans.codingType());
        assertEquals(-2, trans.ExonUpstreamPhase);
        assertEquals(-2, trans.ExonDownstreamPhase);
        assertEquals(8, trans.codingBases()); // stop codon is taken out
        assertEquals(8, trans.totalCodingBases());

        //
        codingStart = new Integer(105);
        codingEnd = new Integer(405);

        transData = createTransExons(
                geneId, transId++, strand, exonStarts, 10, codingStart, codingEnd, true, BIOTYPE_PROTEIN_CODING);

        position = 350;
        trans = extractTranscriptExonData(transData, position, genePosStrand);

        assertEquals(3, trans.ExonUpstream);
        assertEquals(4, trans.ExonDownstream);
        assertEquals(TRANS_CODING_TYPE_CODING, trans.codingType());
        assertEquals(1, trans.ExonUpstreamPhase);
        assertEquals(1, trans.ExonDownstreamPhase);
        assertEquals(28, trans.codingBases());
        assertEquals(31, trans.totalCodingBases());

        // test the reverse strand
        strand = -1;
        geneName = "GENE2";
        geneId = "ENSG0002";
        chromosome = "1";

        GeneAnnotation geneNegStrand = GeneTestUtils.createGeneAnnotation(0, true, geneName, geneId, strand, chromosome, 0, 1);

        // coding taking up exactly the first exon
        codingStart = new Integer(500);
        codingEnd = new Integer(510);

        transData = createTransExons(
                geneId, transId++, strand, exonStarts, 10, codingStart, codingEnd, true, BIOTYPE_PROTEIN_CODING);


        position = 450;
        trans = extractTranscriptExonData(transData, position, geneNegStrand);

        assertEquals(5, trans.ExonMax);
        assertEquals(1, trans.ExonUpstream);
        assertEquals(2, trans.ExonDownstream);
        assertEquals(TRANS_CODING_TYPE_3P_UTR, trans.codingType());
        assertEquals(-2, trans.ExonUpstreamPhase);
        assertEquals(-2, trans.ExonDownstreamPhase);
        assertEquals(8, trans.codingBases());
        assertEquals(8, trans.totalCodingBases());

        //
        codingStart = new Integer(205);
        codingEnd = new Integer(505);

        transData = createTransExons(
                geneId, transId++, strand, exonStarts, 10, codingStart, codingEnd, true, BIOTYPE_PROTEIN_CODING);

        position = 250;
        trans = extractTranscriptExonData(transData, position, geneNegStrand);

        assertEquals(3, trans.ExonUpstream);
        assertEquals(4, trans.ExonDownstream);
        assertEquals(TRANS_CODING_TYPE_CODING, trans.codingType());
        assertEquals(1, trans.ExonUpstreamPhase);
        assertEquals(1, trans.ExonDownstreamPhase);
        assertEquals(28, trans.codingBases());
        assertEquals(31, trans.totalCodingBases());
    }

    @Test
    public void testProteinDomainPositions()
    {
        String geneId = "G0001";
        int transId = 1;

        byte strand = (byte)1;

        int[] exonStarts = new int[]{100, 300, 500};
        int[] exonPhases = new int[]{0, 0, -1};

        TranscriptData transData = createTransExons(geneId, transId++, strand, exonStarts, exonPhases, 100);

        TranscriptProteinData proteinData = new TranscriptProteinData(transId, 0, 0, 5, 55, "hd");

        Integer[] domainPositions = getProteinDomainPositions(proteinData, transData);
        assertEquals(165, (long)domainPositions[SE_START]);
        assertEquals(515, (long)domainPositions[SE_END]);

        // test again with a protein which starts after the first coding exon
        proteinData = new TranscriptProteinData(transId, 0, 0, 55, 65, "hd");

        domainPositions = getProteinDomainPositions(proteinData, transData);
        assertEquals(515, (long)domainPositions[SE_START]);
        assertEquals(545, (long)domainPositions[SE_END]);

        // now on the reverse strand
        strand = (byte)-1;

        proteinData = new TranscriptProteinData(transId, 0, 0, 5, 55, "hd");

        exonStarts = new int[]{100, 300, 500};
        exonPhases = new int[]{-1, 0, 0};

        transData = createTransExons(geneId, transId++, strand, exonStarts, exonPhases, 100);

        domainPositions = getProteinDomainPositions(proteinData, transData);
        assertEquals(185, (long)domainPositions[SE_START]);
        assertEquals(535, (long)domainPositions[SE_END]);

        proteinData = new TranscriptProteinData(transId, 0, 0, 55, 65, "hd");

        domainPositions = getProteinDomainPositions(proteinData, transData);
        assertEquals(155, (long)domainPositions[SE_START]);
        assertEquals(185, (long)domainPositions[SE_END]);

    }

}
