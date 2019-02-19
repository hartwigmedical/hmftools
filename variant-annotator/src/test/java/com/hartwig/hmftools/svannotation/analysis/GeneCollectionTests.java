package com.hartwig.hmftools.svannotation.analysis;

import static com.hartwig.hmftools.svannotation.SvGeneTranscriptCollection.EXON_RANK_MAX;
import static com.hartwig.hmftools.svannotation.SvGeneTranscriptCollection.EXON_RANK_MIN;
import static com.hartwig.hmftools.svannotation.analysis.SvAnnotatorTestUtils.addGeneData;
import static com.hartwig.hmftools.svannotation.analysis.SvAnnotatorTestUtils.addTransExonData;
import static com.hartwig.hmftools.svannotation.analysis.SvAnnotatorTestUtils.createEnsemblGeneData;

import static org.junit.Assert.assertEquals;

import java.util.List;
import java.util.Map;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.variant.structural.annotation.EnsemblGeneData;
import com.hartwig.hmftools.common.variant.structural.annotation.TranscriptExonData;
import com.hartwig.hmftools.svannotation.SvGeneTranscriptCollection;

import org.junit.Test;

public class GeneCollectionTests
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

        int[] transUpExonData = geneTransCache.getExonData(geneName, 11400);

        assertEquals(1, transUpExonData[EXON_RANK_MIN]);
        assertEquals(2, transUpExonData[EXON_RANK_MAX]);

        // before the first
        transUpExonData = geneTransCache.getExonData(geneName, 9000);

        assertEquals(0, transUpExonData[EXON_RANK_MIN]);
        assertEquals(1, transUpExonData[EXON_RANK_MAX]);

        // after the last
        transUpExonData = geneTransCache.getExonData(geneName, 16000);

        assertEquals(4, transUpExonData[EXON_RANK_MIN]);
        assertEquals(-1, transUpExonData[EXON_RANK_MAX]);


        // on an exon boundary
        transUpExonData = geneTransCache.getExonData(geneName, 12500);

        assertEquals(3, transUpExonData[EXON_RANK_MIN]);
        assertEquals(3, transUpExonData[EXON_RANK_MAX]);

    }


}
