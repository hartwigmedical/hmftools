package com.hartwig.hmftools.linx.gene;

import static com.hartwig.hmftools.common.ensemblcache.EnsemblDataCache.EXON_RANK_MAX;
import static com.hartwig.hmftools.common.ensemblcache.EnsemblDataCache.EXON_RANK_MIN;
import static com.hartwig.hmftools.common.ensemblcache.EnsemblDataCache.getProteinDomainPositions;
import static com.hartwig.hmftools.common.gene.GeneTestUtils.createGeneDataCache;
import static com.hartwig.hmftools.common.gene.GeneTestUtils.createTransExons;
import static com.hartwig.hmftools.common.fusion.FusionCommon.NEG_STRAND;
import static com.hartwig.hmftools.common.fusion.FusionCommon.POS_STRAND;
import static com.hartwig.hmftools.common.utils.sv.StartEndIterator.SE_END;
import static com.hartwig.hmftools.common.utils.sv.StartEndIterator.SE_START;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertTrue;

import java.util.List;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.ensemblcache.EnsemblDataCache;
import com.hartwig.hmftools.common.gene.GeneData;
import com.hartwig.hmftools.common.gene.GeneTestUtils;
import com.hartwig.hmftools.common.gene.TranscriptData;
import com.hartwig.hmftools.common.gene.TranscriptProteinData;

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

        List<GeneData> geneList = Lists.newArrayList();
        geneList.add(GeneTestUtils.createEnsemblGeneData(geneId, geneName, chromosome, 1, 10000, 20000));
        GeneTestUtils.addGeneData(geneTransCache, chromosome, geneList);

        List<TranscriptData> transDataList = Lists.newArrayList();

        int transId = 1;
        byte strand = 1;

        int[] exonStarts = new int[]{10500, 11500, 12500, 13500};

        int codingStart = 11500;
        int codingEnd = 13598;
        TranscriptData transData = createTransExons(geneId, transId, strand, exonStarts, 99, codingStart, codingEnd,true, "");
        transDataList.add(transData);

        int transId2 = 2;

        exonStarts = new int[]{12500, 13500, 14500};

        transData = createTransExons(geneId, transId2, strand, exonStarts, 100, null, null, true, "");
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
    public void testProteinDomainPositions()
    {
        String geneId = "G0001";
        int transId = 1;

        int[] exonStarts = new int[]{100, 300, 500};

        int codingStart = 150;
        int codingEnd = 550;
        TranscriptData transData = createTransExons(geneId, transId++, POS_STRAND, exonStarts, 100, codingStart, codingEnd, true, "");

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
        proteinData = new TranscriptProteinData(transId, 0, 0, 5, 55, "hd");

        exonStarts = new int[]{100, 300, 500};

        codingStart = 350;
        codingEnd = 550;
        transData = createTransExons(geneId, transId++, NEG_STRAND, exonStarts, 100, codingStart, codingEnd, true, "");

        domainPositions = getProteinDomainPositions(proteinData, transData);
        assertEquals(185, (long)domainPositions[SE_START]);
        assertEquals(535, (long)domainPositions[SE_END]);

        proteinData = new TranscriptProteinData(transId, 0, 0, 55, 65, "hd");

        domainPositions = getProteinDomainPositions(proteinData, transData);
        assertEquals(155, (long)domainPositions[SE_START]);
        assertEquals(185, (long)domainPositions[SE_END]);

    }

}
