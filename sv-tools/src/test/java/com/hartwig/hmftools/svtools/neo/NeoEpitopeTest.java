package com.hartwig.hmftools.svtools.neo;

import static com.hartwig.hmftools.linx.fusion.FusionConstants.PRE_GENE_PROMOTOR_DISTANCE;
import static com.hartwig.hmftools.common.ensemblcache.GeneTestUtils.createGeneDataCache;
import static com.hartwig.hmftools.common.ensemblcache.GeneTestUtils.addGeneData;
import static com.hartwig.hmftools.common.ensemblcache.GeneTestUtils.addTransExonData;
import static com.hartwig.hmftools.common.ensemblcache.GeneTestUtils.createEnsemblGeneData;
import static com.hartwig.hmftools.common.ensemblcache.GeneTestUtils.createTransExons;
import static com.hartwig.hmftools.common.ensemblcache.GeneTestUtils.generateTransName;
import static com.hartwig.hmftools.svtools.neo.AminoAcidConverter.reverseStrandBases;
import static com.hartwig.hmftools.svtools.neo.AminoAcidConverter.swapDnaToRna;
import static com.hartwig.hmftools.svtools.neo.AminoAcidConverter.swapRnaToDna;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertTrue;

import java.util.List;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.ensemblcache.EnsemblDataCache;
import com.hartwig.hmftools.common.ensemblcache.EnsemblGeneData;
import com.hartwig.hmftools.common.ensemblcache.ExonData;
import com.hartwig.hmftools.common.fusion.GeneAnnotation;
import com.hartwig.hmftools.common.ensemblcache.TranscriptData;
import com.hartwig.hmftools.linx.fusion.FusionParameters;
import com.hartwig.hmftools.linx.fusion.GeneFusion;
import com.hartwig.hmftools.linx.fusion.NeoEpitopeWriter;
import com.hartwig.hmftools.common.genome.refgenome.MockRefGenome;

import org.junit.Test;

public class NeoEpitopeTest
{
    @Test
    public void testDnaRnaRoutines()
    {
        String dnaBases = "AGCT";
        String rnaBases = swapDnaToRna(dnaBases);
        assertTrue(rnaBases.equals("AGCU"));
        assertTrue(dnaBases.equals(swapRnaToDna(rnaBases)));

        dnaBases = "AGCTTCGACT";
        String reverseStrandDna = reverseStrandBases(dnaBases);
        assertTrue(reverseStrandDna.equals("AGTCGAAGCT"));
    }

    /*
    @Test
    public void testNeoEpitopes()
    {
        LinxTester tester = new LinxTester();

        EnsemblDataCache geneTransCache = createGeneDataCache();

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

        int[] exonStarts = new int[]{5, 15, 25, 35, 45, 55, 65, 75};
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

        NeoEpitopeWriter neoEpFinder = new NeoEpitopeWriter(refGenome, geneTransCache, "");

        byte posOrient = 1;
        byte negOrient = -1;

        // create fusions between various phases sections of these transcripts
        // add upstream breakends
        List<GeneAnnotation> upGenes = Lists.newArrayList();
        int upPos = 22;
        upGenes.addAll(geneTransCache.findGeneAnnotationsBySv(0, true, chromosome1, upPos, posOrient, PRE_GENE_PROMOTOR_DISTANCE));
        upGenes.get(0).setPositionalData(chromosome1, upPos, posOrient);

        // add downstream breakends
        List<GeneAnnotation> downGenes = Lists.newArrayList();
        int downPos = 22;
        downGenes.addAll(geneTransCache.findGeneAnnotationsBySv(0, false, chromosome2, downPos, negOrient, PRE_GENE_PROMOTOR_DISTANCE));
        downGenes.get(0).setPositionalData(chromosome2, downPos, negOrient);

        FusionParameters params = new FusionParameters();

        // List<GeneFusion> fusions = tester.FusionAnalyser.getFusionFinder().findFusions(upGenes, downGenes, params, false);
        List<GeneFusion> fusions = tester.FusionAnalyser.getFusionFinder().findFusions(upGenes, downGenes, params);

        neoEpFinder.reportNeoEpitopes(tester.SampleId, fusions);

        assertEquals(1, neoEpFinder.getResults().size());
        NeoEpitopeFusion data = neoEpFinder.getResults().get(0);
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

        fusions = tester.FusionAnalyser.getFusionFinder().findFusions(upGenes, downGenes,params);

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

        fusions = tester.FusionAnalyser.getFusionFinder().findFusions(upGenes, downGenes, params);

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

        fusions = tester.FusionAnalyser.getFusionFinder().findFusions(upGenes, downGenes, params);

        neoEpFinder.reportNeoEpitopes(tester.SampleId, fusions);

        assertEquals(1, neoEpFinder.getResults().size());
        data = neoEpFinder.getResults().get(0);
        assertTrue(data.upstreamAcids().equals("MSSS"));
        assertTrue(data.novelAcid().equals("YHHHHHHH"));
        assertTrue(data.downstreamAcids().equals(""));
    }

    @Test
    public void testNeoEpitopesReverseStrand()
    {
        LinxTester tester = new LinxTester();

        EnsemblDataCache geneTransCache = createGeneDataCache();

        tester.initialiseFusions(geneTransCache);

        PRE_GENE_PROMOTOR_DISTANCE = 10;

        // first a gene on the forward strand
        String geneId1 = "ENSG0001";
        String chromosome1 = "1";
        byte strand = -1;

        int transId = 1;

        String intron = "AAAAA";
        String nonCodingExon = "GGGGG";
        String nonGenicDna = "TTTTT";
        String sCodon = convertAminoAcidToDnaCodon("S");

        // resultant bases and exon indices
        // AAAAACCCCCTTATGTTTTTATGATGATGATGATGATGATGATGATGATGATTTTTTGACATCCCCCAAAAA
        // 0    5    10   15   20                             51   56    62   67

        // first exon, starts coding
        String exon1 = nonCodingExon + swapRnaToDna(START_CODON) + sCodon;
        String transBases = exon1;

        // second exon - a string of 'S' amino acids following by 1 single base (ending on phase 1)
        String exon2 = "";
        int codonCount = 10;
        for(int i = 0; i < codonCount; ++i)
        {
            exon2 += sCodon;
        }

        // finish on phase 1
        exon2 += sCodon.substring(0, 1);

        transBases += intron + exon2;

        // final exon
        String exon3 = sCodon.substring(1) + swapRnaToDna(STOP_CODON_1) + nonCodingExon;
        transBases += intron + exon3;

        String revBases = reverseStrandBases(nonGenicDna + transBases + nonGenicDna);

        int transStart = nonGenicDna.length();
        int transEnd = transStart + transBases.length() - 1;

        // exons added from lower to higher positions
        List<ExonData> exons = Lists.newArrayList();
        exons.add(new ExonData(transId, transStart, transStart + exon3.length() - 1, 3, 1, -1));

        int nextExonStart = exons.get(exons.size() - 1).ExonEnd + intron.length() + 1;
        exons.add(new ExonData(transId, nextExonStart, nextExonStart + exon2.length() - 1, 2, 0, 1));

        nextExonStart = exons.get(exons.size() - 1).ExonEnd + intron.length() + 1;
        exons.add(new ExonData(transId, nextExonStart, nextExonStart + exon1.length() - 1, 1, -1, 0));

        TranscriptData transData = new TranscriptData(transId, generateTransName(transId), geneId1, true, strand, transStart, transEnd,
                transStart + nonCodingExon.length(), transEnd - nonCodingExon.length(), "");

        transData.exons().addAll(exons);

        addTransExonData(geneTransCache, geneId1, Lists.newArrayList(transData));

        // same format transcript on another chromosome, same strand
        String geneId2 = "ENSG0002";
        String chromosome2 = "2";

        List<EnsemblGeneData> geneList = Lists.newArrayList();
        geneList.add(createEnsemblGeneData(geneId1, "GENE1", chromosome1, strand, transStart, transEnd));
        addGeneData(geneTransCache, chromosome1, geneList);

        geneList = Lists.newArrayList();
        geneList.add(createEnsemblGeneData(geneId2, "GENE2", chromosome2, strand, transStart, transEnd));
        addGeneData(geneTransCache, chromosome2, geneList);

        transData = new TranscriptData(++transId, generateTransName(transId), geneId2, true, strand, transStart, transEnd,
                transStart + nonCodingExon.length(), transEnd - nonCodingExon.length(), "");

        transData.exons().addAll(exons);

        addTransExonData(geneTransCache, geneId2, Lists.newArrayList(transData));

        MockRefGenome refGenome = new MockRefGenome();

        refGenome.RefGenomeMap.put(chromosome1, revBases);
        refGenome.RefGenomeMap.put(chromosome2, revBases);

        NeoEpitopeWriter neoEpFinder = new NeoEpitopeWriter(refGenome, geneTransCache, "");

        byte posOrient = 1;
        byte negOrient = -1;

        // create fusions between various phases sections of these transcripts
        // add upstream breakends
        List<GeneAnnotation> upGenes = Lists.newArrayList();
        int upPos = 51;
        upGenes.addAll(geneTransCache.findGeneAnnotationsBySv(0, true, chromosome1, upPos, posOrient, PRE_GENE_PROMOTOR_DISTANCE));
        upGenes.get(0).setPositionalData(chromosome1, upPos, posOrient);

        // add downstream breakends
        List<GeneAnnotation> downGenes = Lists.newArrayList();
        int downPos = 51;
        downGenes.addAll(geneTransCache.findGeneAnnotationsBySv(0, false, chromosome2, downPos, negOrient, PRE_GENE_PROMOTOR_DISTANCE));
        downGenes.get(0).setPositionalData(chromosome2, downPos, negOrient);

        FusionParameters params = new FusionParameters();

        List<GeneFusion> fusions = tester.FusionAnalyser.getFusionFinder().findFusions(upGenes, downGenes, params);

        neoEpFinder.reportNeoEpitopes(tester.SampleId, fusions);

        assertEquals(1, neoEpFinder.getResults().size());
        NeoEpitopeFusion data = neoEpFinder.getResults().get(0);
        assertTrue(data.upstreamAcids().equals("MS"));
        assertTrue(data.novelAcid().equals(""));
        assertTrue(data.downstreamAcids().equals("SSSSSSSSSS"));

        // with a fusion between the 2nd and 3rd exons, splitting a codon
        upPos = 17;
        upGenes.clear();
        upGenes.addAll(geneTransCache.findGeneAnnotationsBySv(0, true, chromosome1, upPos, posOrient, PRE_GENE_PROMOTOR_DISTANCE));
        upGenes.get(0).setPositionalData(chromosome1, upPos, posOrient);

        // add downstream breakends
        downGenes.clear();
        downPos = 17;
        downGenes.addAll(geneTransCache.findGeneAnnotationsBySv(0, false, chromosome2, downPos, negOrient, PRE_GENE_PROMOTOR_DISTANCE));
        downGenes.get(0).setPositionalData(chromosome2, downPos, negOrient);

        fusions = tester.FusionAnalyser.getFusionFinder().findFusions(upGenes, downGenes, params);

        neoEpFinder.reportNeoEpitopes(tester.SampleId, fusions);

        assertEquals(1, neoEpFinder.getResults().size());
        data = neoEpFinder.getResults().get(0);
        assertTrue(data.upstreamAcids().equals("SSSSSSSSSS"));
        assertTrue(data.novelAcid().equals("S"));
        assertTrue(data.downstreamAcids().equals("_"));

        // a fusion to the 5'UTR of the downstream gene
        upPos = 51;
        upGenes.clear();
        upGenes.addAll(geneTransCache.findGeneAnnotationsBySv(0, true, chromosome1, upPos, negOrient, PRE_GENE_PROMOTOR_DISTANCE));
        upGenes.get(0).setPositionalData(chromosome1, upPos, negOrient);

        // add downstream breakends
        downGenes.clear();
        downPos = 68;
        downGenes.addAll(geneTransCache.findGeneAnnotationsBySv(0, false, chromosome2, downPos, posOrient, PRE_GENE_PROMOTOR_DISTANCE));
        downGenes.get(0).setPositionalData(chromosome2, downPos, posOrient);

        fusions = tester.FusionAnalyser.getFusionFinder().findFusions(upGenes, downGenes, params);

        neoEpFinder.reportNeoEpitopes(tester.SampleId, fusions);

        assertEquals(1, neoEpFinder.getResults().size());
        data = neoEpFinder.getResults().get(0);
        assertTrue(data.upstreamAcids().equals("MS"));
        assertTrue(data.novelAcid().equals(""));
        assertTrue(data.downstreamAcids().equals("SSSSSSSSSS"));

        // a fusion to the 3'UTR of the downstream gene
        upPos = 51;
        upGenes.clear();
        upGenes.addAll(geneTransCache.findGeneAnnotationsBySv(0, true, chromosome1, upPos, negOrient, PRE_GENE_PROMOTOR_DISTANCE));
        upGenes.get(0).setPositionalData(chromosome1, upPos, negOrient);

        // add downstream breakends
        downGenes.clear();
        downPos = 14;
        downGenes.addAll(geneTransCache.findGeneAnnotationsBySv(0, false, chromosome2, downPos, posOrient, PRE_GENE_PROMOTOR_DISTANCE));
        downGenes.get(0).setPositionalData(chromosome2, downPos, posOrient);

        fusions.clear();

        neoEpFinder.checkFusions(fusions, upGenes, downGenes);
        assertEquals(1, fusions.size());

        neoEpFinder.reportNeoEpitopes(tester.SampleId, fusions);

        assertEquals(1, neoEpFinder.getResults().size());
        data = neoEpFinder.getResults().get(0);
        assertTrue(data.upstreamAcids().equals("MS"));
        assertTrue(data.novelAcid().equals("H"));
        assertTrue(data.downstreamAcids().equals(""));

        // again with a longer downstream 3'UTR sequence
        String geneId3 = "ENSG0003";
        String chromosome3 = "3";
        strand = 1;

        // resultant bases and exon indices
        // AAAAACCCCCTTATGTTTTTATGATGATGATGATGATGATGATGATGATGATTTTTTGACATCCCCCAAAAA
        // 0    5    10   15   20                             51   56    62   67

        transBases = nonCodingExon; // exon 1
        transBases += intron + nonCodingExon + swapRnaToDna(START_CODON) + sCodon; // exon 2
        transBases += intron + sCodon + swapRnaToDna(STOP_CODON_1) + nonCodingExon; // exon 3
        transBases += intron + nonCodingExon; // exon 4
        transBases += intron + nonCodingExon; // exon 5
        transBases += intron + nonCodingExon; // exon 6

        String refBases = nonGenicDna + transBases + nonGenicDna;

        refGenome.RefGenomeMap.put(chromosome3, refBases);

        transStart = nonGenicDna.length();
        transEnd = transStart + transBases.length() - 1;


        exons.clear();
        exons.add(new ExonData(transId, transStart, transStart + nonCodingExon.length() - 1, 1, -1, -1));

        nextExonStart = exons.get(exons.size() - 1).ExonEnd + intron.length() + 1;
        exons.add(new ExonData(transId, nextExonStart, nextExonStart + nonCodingExon.length() + 6 - 1, 2, -1, 0));

        nextExonStart = exons.get(exons.size() - 1).ExonEnd + intron.length() + 1;
        exons.add(new ExonData(transId, nextExonStart, nextExonStart + 6 + nonCodingExon.length() - 1, 3, 0, -1));

        nextExonStart = exons.get(exons.size() - 1).ExonEnd + intron.length() + 1;
        exons.add(new ExonData(transId, nextExonStart, nextExonStart + nonCodingExon.length() - 1, 4, -1, -1));

        nextExonStart = exons.get(exons.size() - 1).ExonEnd + intron.length() + 1;
        exons.add(new ExonData(transId, nextExonStart, nextExonStart + nonCodingExon.length() - 1, 5, -1, -1));

        nextExonStart = exons.get(exons.size() - 1).ExonEnd + intron.length() + 1;
        exons.add(new ExonData(transId, nextExonStart, nextExonStart + nonCodingExon.length() - 1, 6, -1, -1));

        int codingStart = exons.get(1).ExonStart + nonCodingExon.length() - 1;
        int codingEnd = exons.get(2).ExonStart + 6 - 1;

        transData = new TranscriptData(++transId, generateTransName(transId), geneId3, true, strand, transStart, transEnd,
                codingStart, codingEnd, "");

        transData.exons().addAll(exons);

        geneList = Lists.newArrayList();
        geneList.add(createEnsemblGeneData(geneId3, "GENE3", chromosome3, strand, transStart, transEnd));
        addGeneData(geneTransCache, chromosome3, geneList);

        addTransExonData(geneTransCache, geneId3, Lists.newArrayList(transData));

        // a fusion to the 3'UTR of the downstream gene
        upPos = 51;
        upGenes.clear();
        upGenes.addAll(geneTransCache.findGeneAnnotationsBySv(0, true, chromosome1, upPos, negOrient, PRE_GENE_PROMOTOR_DISTANCE));
        upGenes.get(0).setPositionalData(chromosome1, upPos, negOrient);

        // add downstream breakends
        downGenes.clear();
        downPos = codingEnd + nonCodingExon.length() + 2;
        downGenes.addAll(geneTransCache.findGeneAnnotationsBySv(0, false, chromosome3, downPos, negOrient, PRE_GENE_PROMOTOR_DISTANCE));
        downGenes.get(0).setPositionalData(chromosome3, downPos, negOrient);

        fusions.clear();

        neoEpFinder.checkFusions(fusions, upGenes, downGenes);
        assertEquals(1, fusions.size());

        neoEpFinder.reportNeoEpitopes(tester.SampleId, fusions);

        assertEquals(1, neoEpFinder.getResults().size());
        data = neoEpFinder.getResults().get(0);
        assertTrue(data.upstreamAcids().equals("MS"));
        assertTrue(data.novelAcid().equals("GGGGG"));
        assertTrue(data.downstreamAcids().equals(""));
        assertEquals(10, data.downstreamNmdBases());
    }
    */

}
