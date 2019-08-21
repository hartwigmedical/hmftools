package com.hartwig.hmftools.linx.fusion;

import static com.hartwig.hmftools.linx.gene.SvGeneTranscriptCollection.PRE_GENE_PROMOTOR_DISTANCE;
import static com.hartwig.hmftools.linx.utils.GeneTestUtils.addGeneData;
import static com.hartwig.hmftools.linx.utils.GeneTestUtils.addTransExonData;
import static com.hartwig.hmftools.linx.utils.GeneTestUtils.createEnsemblGeneData;
import static com.hartwig.hmftools.linx.utils.GeneTestUtils.createTransExons;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertTrue;

import java.util.List;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.variant.structural.annotation.EnsemblGeneData;
import com.hartwig.hmftools.common.variant.structural.annotation.GeneAnnotation;
import com.hartwig.hmftools.common.variant.structural.annotation.GeneFusion;
import com.hartwig.hmftools.common.variant.structural.annotation.TranscriptData;
import com.hartwig.hmftools.linx.gene.SvGeneTranscriptCollection;
import com.hartwig.hmftools.linx.neoepitope.NeoEpitopeData;
import com.hartwig.hmftools.linx.neoepitope.NeoEpitopeFinder;
import com.hartwig.hmftools.linx.neoepitope.RefGenomeSource;
import com.hartwig.hmftools.linx.utils.LinxTester;
import com.hartwig.hmftools.linx.utils.MockRefGenome;

import org.junit.Test;

public class NeoEpitopeTest
{
    @Test
    public void testNeoEpitopes()
    {
        LinxTester tester = new LinxTester();
        tester.logVerbose(true);

        SvGeneTranscriptCollection geneTransCache = new SvGeneTranscriptCollection();

        tester.initialiseFusions(geneTransCache);

        PRE_GENE_PROMOTOR_DISTANCE = 10;

        // first a gene on the forward strand
        String geneName = "GENE1";
        String geneId1 = "ENSG0001";
        String chromosome1 = "1";
        byte strand = 1;

        List<EnsemblGeneData> geneList = Lists.newArrayList();
        geneList.add(createEnsemblGeneData(geneId1, geneName, chromosome1, strand, 5, 100));
        addGeneData(geneTransCache, chromosome1, geneList);

        int transId = 1;

        long[] exonStarts = new long[]{5, 15, 25, 35, 45, 55, 65, 75};
        int[] exonPhases = new int[]{-1, 0, 2, 1, 0, 2, 1, -1};

        TranscriptData transData = createTransExons(geneId1, transId++, strand, exonStarts, exonPhases, 4, true);
        assertEquals(17, transData.CodingStart.longValue());

        List<TranscriptData> transDataList = Lists.newArrayList(transData);

        addTransExonData(geneTransCache, geneId1, transDataList);

        geneName = "GENE2";
        String geneId2 = "ENSG0002";
        String chromosome2 = "2";

        geneList = Lists.newArrayList();
        geneList.add(createEnsemblGeneData(geneId2, geneName, chromosome2, strand, 5, 100));

        addGeneData(geneTransCache, chromosome2, geneList);

        transData = createTransExons(geneId2, transId++, strand, exonStarts, exonPhases, 4, true);
        transDataList = Lists.newArrayList(transData);

        addTransExonData(geneTransCache, geneId2, transDataList);

        MockRefGenome refGenome = new MockRefGenome();

        String refBases = "";
        String intron = "AAAAA";
        String nonCodingExon = "GGGGG";

        refBases += intron;
        refBases += nonCodingExon; // exon 1
        refBases += intron;
        refBases += "GGATG"; // exon 2, start of coding, ends on phase 2, so next is 0
        refBases += intron;
        refBases += "TCATC"; // exon 3, end phase 1
        refBases += intron;
        refBases += "ATCAT"; // exon 4, end phase 0
        refBases += intron;
        refBases += "CATCA"; // exon 5, end phase 2
        refBases += intron;
        refBases += "TCATC"; // exon 6, end phase 1
        refBases += intron;
        refBases += "ATCAT"; // exon 7, end phase 0
        refBases += intron;
        refBases += "CATAA"; // exon 8 including stop codon
        refBases += intron + intron;

        refGenome.RefGenomeMap.put(chromosome1, refBases);
        refGenome.RefGenomeMap.put(chromosome2, refBases);

        NeoEpitopeFinder neoEpFinder = new NeoEpitopeFinder(refGenome, geneTransCache, "");

        byte posOrient = 1;
        byte negOrient = -1;

        // create fusions between various phases sections of these transcripts
        // add upstream breakends
        List<GeneAnnotation> upGenes = Lists.newArrayList();
        long upPos = 22;
        upGenes.addAll(geneTransCache.findGeneAnnotationsBySv(0, true, chromosome1, upPos, posOrient, PRE_GENE_PROMOTOR_DISTANCE));
        upGenes.get(0).setPositionalData(chromosome1, upPos, posOrient);

        // add downstream breakends
        List<GeneAnnotation> downGenes = Lists.newArrayList();
        long downPos = 22;
        downGenes.addAll(geneTransCache.findGeneAnnotationsBySv(0, false, chromosome2, downPos, negOrient, PRE_GENE_PROMOTOR_DISTANCE));
        downGenes.get(0).setPositionalData(chromosome2, downPos, negOrient);

        List<GeneFusion> fusions = tester.FusionAnalyser.getFusionFinder().findFusions(upGenes, downGenes,
                false, true, null, false);

        neoEpFinder.reportNeoEpitopes(tester.SampleId, fusions);

        assertEquals(1, neoEpFinder.getResults().size());
        NeoEpitopeData data = neoEpFinder.getResults().get(0);
        assertTrue(data.upstreamAcids().equals("M"));
        assertTrue(data.novelAcid().equals(""));
        assertTrue(data.downstreamAcids().equals("SSSSSSSSS"));

        // try again a phase 1 fusion
        upPos = 41;
        upGenes.clear();
        upGenes.addAll(geneTransCache.findGeneAnnotationsBySv(0, true, chromosome1, upPos, posOrient, PRE_GENE_PROMOTOR_DISTANCE));
        upGenes.get(0).setPositionalData(chromosome1, upPos, posOrient);

        // add downstream breakends
        downPos = 41;
        downGenes.clear();
        downGenes.addAll(geneTransCache.findGeneAnnotationsBySv(0, false, chromosome2, downPos, negOrient, PRE_GENE_PROMOTOR_DISTANCE));
        downGenes.get(0).setPositionalData(chromosome2, downPos, negOrient);

        fusions = tester.FusionAnalyser.getFusionFinder().findFusions(upGenes, downGenes,
                false, true, null, false);

        neoEpFinder.reportNeoEpitopes(tester.SampleId, fusions);

        assertEquals(1, neoEpFinder.getResults().size());
        data = neoEpFinder.getResults().get(0);
        assertTrue(data.upstreamAcids().equals("MSSS"));
        assertTrue(data.novelAcid().equals("S"));
        assertTrue(data.downstreamAcids().equals("SSSSS"));

        // phase 2 fusion
        upPos = 62;
        upGenes.clear();
        upGenes.addAll(geneTransCache.findGeneAnnotationsBySv(0, true, chromosome1, upPos, posOrient, PRE_GENE_PROMOTOR_DISTANCE));
        upGenes.get(0).setPositionalData(chromosome1, upPos, posOrient);

        // add downstream breakends
        downPos = 62;
        downGenes.clear();
        downGenes.addAll(geneTransCache.findGeneAnnotationsBySv(0, false, chromosome2, downPos, negOrient, PRE_GENE_PROMOTOR_DISTANCE));
        downGenes.get(0).setPositionalData(chromosome2, downPos, negOrient);

        fusions = tester.FusionAnalyser.getFusionFinder().findFusions(upGenes, downGenes,
                false, true, null, false);

        neoEpFinder.reportNeoEpitopes(tester.SampleId, fusions);

        assertEquals(1, neoEpFinder.getResults().size());
        data = neoEpFinder.getResults().get(0);
        assertTrue(data.upstreamAcids().equals("MSSSSSS"));
        assertTrue(data.novelAcid().equals("S"));
        assertTrue(data.downstreamAcids().equals("SS"));

        // unphased fusion - only difference is that all downstream bases are collected (and may not form viable amino acids)
        upPos = 41;
        upGenes.clear();
        upGenes.addAll(geneTransCache.findGeneAnnotationsBySv(0, true, chromosome1, upPos, posOrient, PRE_GENE_PROMOTOR_DISTANCE));
        upGenes.get(0).setPositionalData(chromosome1, upPos, posOrient);

        // add downstream breakends
        downPos = 31;
        downGenes.clear();
        downGenes.addAll(geneTransCache.findGeneAnnotationsBySv(0, false, chromosome2, downPos, negOrient, PRE_GENE_PROMOTOR_DISTANCE));
        downGenes.get(0).setPositionalData(chromosome2, downPos, negOrient);

        fusions = tester.FusionAnalyser.getFusionFinder().findFusions(upGenes, downGenes,
                false, true, null, false);

        neoEpFinder.reportNeoEpitopes(tester.SampleId, fusions);

        assertEquals(1, neoEpFinder.getResults().size());
        data = neoEpFinder.getResults().get(0);
        assertTrue(data.upstreamAcids().equals("MSSS"));
        assertTrue(data.novelAcid().equals("Y"));
        assertTrue(data.downstreamAcids().equals("HHHHHHH"));
    }

}
