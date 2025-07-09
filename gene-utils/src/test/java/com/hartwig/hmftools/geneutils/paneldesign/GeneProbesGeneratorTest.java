package com.hartwig.hmftools.geneutils.paneldesign;

import static com.hartwig.hmftools.common.test.GeneTestUtils.CHR_1;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertTrue;
import static org.junit.Assert.fail;

import java.util.Iterator;
import java.util.List;

import com.hartwig.hmftools.common.gene.ExonData;
import com.hartwig.hmftools.common.gene.GeneData;
import com.hartwig.hmftools.common.gene.TranscriptData;
import com.hartwig.hmftools.common.genome.region.Strand;

import com.hartwig.hmftools.common.region.BaseRegion;
import com.hartwig.hmftools.common.test.MockRefGenome;

import org.apache.commons.lang3.StringUtils;
import org.junit.Test;

public class GeneProbesGeneratorTest
{
    @Test
    public void testPopulateTargetGeneRegionsSimple()
    {
//        fail();
//        // create gene and transcript data
//        String geneId = "GENE01";
//
//        GeneData geneData = new GeneData(geneId, geneId, CHR_1, Strand.POS_STRAND, 11000, 19000 - 1, "");
//
//        int transId = 1;
//        String transName = "TRANS01";
//
//        TranscriptData transData = new TranscriptData(
//                transId, transName, geneId, true, Strand.POS_STRAND,
//                11000, 19000 - 1, 13500, 18000 - 1, "", null);
//
//        // first exon is not coding
//        transData.exons().add(new ExonData(transId, 11000, 12000 - 1, 1, -1, -1));
//
//        // coding exon
//        transData.exons().add(new ExonData(transId, 13000, 14000 - 1, 2, -1, -1));
//
//        // coding exon
//        transData.exons().add(new ExonData(transId, 15000, 19000 - 1, 3, -1, -1));
//
//        // create a targeted gene
//        Gene targetedGene = new Gene(geneData, transData);
//
//        GeneProbesGenerator.populateTargetedGeneRegions(targetedGene);
//
//        // check the regions
//        assertEquals(5, targetedGene.getRegions().size());
//
//        Iterator<GeneRegion> itr = targetedGene.getRegions().iterator();
//
//        // first we should have an upstream region
//        GeneRegion region = itr.next();
//        assertEquals(GeneRegion.Type.UP_STREAM, region.getType());
//        assertEquals(9040, region.getStart());
//        assertEquals(9999, region.getEnd());
//
//        // 1st exon is non coding, add a probe in the centre
//        region = itr.next();
//        assertEquals(GeneRegion.Type.UTR, region.getType());
//        assertEquals(11440, region.getStart());
//        assertEquals(11559, region.getEnd());
//
//        // 2nd exon is coding, cover the coding regions of the exon
//        region = itr.next();
//        assertEquals(GeneRegion.Type.CODING, region.getType());
//        assertEquals(13500, region.getStart());
//        assertEquals(14000 - 1, region.getEnd());
//
//        // 3rd exon is also coding, cover the coding regions of the exon
//        region = itr.next();
//        assertEquals(GeneRegion.Type.CODING, region.getType());
//        assertEquals(15000, region.getStart());
//        assertEquals(18000 - 1, region.getEnd());
//
//        // downstream region
//        region = itr.next();
//        assertEquals(GeneRegion.Type.DOWN_STREAM, region.getType());
//        assertEquals(20000, region.getStart());
//        assertEquals(20960 - 1, region.getEnd());
    }

    @Test
    public void testPopulateTargetGeneRegionsIntrons()
    {
//        fail();
//        // create gene and transcript data
//        String geneId = "GENE01";
//
//        GeneData geneData = new GeneData(
//                geneId, geneId, CHR_1, Strand.POS_STRAND, 11000, 25000 - 1, "");
//
//        int transId = 1;
//        String transName = "TRANS01";
//
//        TranscriptData transData = new TranscriptData(
//                transId, transName, geneId, true, Strand.POS_STRAND,
//                11000, 25000 - 1, 13500, 24000 - 1, "", null);
//
//        // first exon is not coding
//        transData.exons().add(new ExonData(transId, 11000, 12000 - 1, 1, -1, -1));
//
//        // long intron 6kb (>5kb)
//
//        // coding exon
//        transData.exons().add(new ExonData(transId, 18000, 19000 - 1, 2, -1, -1));
//
//        // short intron 4kb (3-5kb)
//
//        // coding exon
//        transData.exons().add(new ExonData(transId, 23000, 25000 - 1, 3, -1, -1));
//
//        // create a targeted gene
//        Gene targetedGene = new Gene(geneData, transData);
//
//        GeneProbesGenerator.populateTargetedGeneRegions(targetedGene);
//
//        // check the regions
//        assertEquals(8, targetedGene.getRegions().size());
//        Iterator<GeneRegion> itr = targetedGene.getRegions().iterator();
//
//        // first we should have an upstream region
//        GeneRegion region = itr.next();
//        assertEquals(GeneRegion.Type.UP_STREAM, region.getType());
//        assertEquals(9040, region.getStart());
//        assertEquals(9999, region.getEnd());
//
//        // 1st exon is non coding, add a probe in the centre
//        region = itr.next();
//        assertEquals(GeneRegion.Type.UTR, region.getType());
//        assertEquals(11440, region.getStart());
//        assertEquals(11559, region.getEnd());
//
//        // long intron create two flanking regions
//        region = itr.next();
//        assertEquals(GeneRegion.Type.INTRONIC_LONG, region.getType());
//        assertEquals(13000, region.getStart()); // 1000 bases after previous exon
//        assertEquals(13960 - 1, region.getEnd());
//
//        region = itr.next();
//        assertEquals(GeneRegion.Type.INTRONIC_LONG, region.getType());
//        assertEquals(16040, region.getStart());
//        assertEquals(17000 - 1, region.getEnd()); // 1000 bases before next exon
//
//        // 2nd exon is coding, cover the coding regions of the exon
//        region = itr.next();
//        assertEquals(GeneRegion.Type.CODING, region.getType());
//        assertEquals(18000, region.getStart());
//        assertEquals(19000 - 1, region.getEnd());
//
//        // short intron, 960 bases length at 21000, which is the middle
//        region = itr.next();
//        assertEquals(GeneRegion.Type.INTRONIC_SHORT, region.getType());
//        assertEquals(20520, region.getStart());
//        assertEquals(21480 - 1, region.getEnd());
//
//        // 3rd exon
//        region = itr.next();
//        assertEquals(GeneRegion.Type.CODING, region.getType());
//        assertEquals(23000, region.getStart());
//        assertEquals(24000 - 1, region.getEnd());
//
//        // downstream region
//        region = itr.next();
//        assertEquals(GeneRegion.Type.DOWN_STREAM, region.getType());
//        assertEquals(26000, region.getStart());
//        assertEquals(26960 - 1, region.getEnd());
    }

    @Test
    public void testPopulateCandidateProbes()
    {
//        fail();
//        MockRefGenome refGenome = new MockRefGenome();
//
//        // make it such that probe GC is 33.333%
//        refGenome.RefGenomeMap.put("X", StringUtils.repeat("ATGTTCAAGTAC", 1000));
//
//        // create gene and transcript data
//        String geneId = "GENE01";
//
//        GeneData geneData = new GeneData(geneId, geneId, "X", Strand.POS_STRAND, 1, 1000, "");
//
//        int transId = 1;
//        String transName = "TRANS01";
//
//        TranscriptData transData = new TranscriptData(
//                transId, transName, geneId, true, Strand.POS_STRAND, 1, 1000, 1, 1000, "", null);
//
//        // create a targeted gene
//        Gene targetedGene = new Gene(geneData, transData);
//
//        // coding gene region
//        GeneRegion region = new GeneRegion(targetedGene, GeneRegion.Type.CODING, new BaseRegion(1, 1000));
//        GeneProbesGenerator.populateCandidateProbes(region, refGenome);
//        // for coding exons, we use the whole region, so no probe candidate is needed
//        assertTrue(region.useWholeRegion());
//        assertTrue(region.getProbeCandidates().isEmpty());
//
//        // test that we generate the correct probe candidates and the gc
//        region = new GeneRegion(targetedGene, GeneRegion.Type.UP_STREAM, new BaseRegion(1, 1000));
//        GeneProbesGenerator.populateCandidateProbes(region, refGenome);
//
//        // 8 probes generated
//        assertEquals(8, region.getProbeCandidates().size());
//
//        checkProbeCandidates(region.getProbeCandidates(), region.getStart(), region.getEnd(), 8, 120);
    }

//    private void checkProbeCandidates(List<CandidateProbe> probeCandidates, int start, int end, int numProbes, int probeLength)
//    {
//        int probeStart = start;
//        for(CandidateProbe probeCandidate : probeCandidates)
//        {
//            assertEquals(probeStart, probeCandidate.getStart());
//            assertEquals(probeStart + probeLength - 1, probeCandidate.getEnd());
//            assertEquals(1.0 / 3.0, probeCandidate.getGcContent(), 1e-5);
//            probeStart += probeLength;
//        }
//    }
}