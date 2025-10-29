package com.hartwig.hmftools.linx.fusion;

import static com.hartwig.hmftools.common.test.GeneTestUtils.addGeneData;
import static com.hartwig.hmftools.common.test.GeneTestUtils.addTransExonData;
import static com.hartwig.hmftools.common.test.GeneTestUtils.createEnsemblGeneData;
import static com.hartwig.hmftools.common.test.GeneTestUtils.createGeneDataCache;
import static com.hartwig.hmftools.common.test.GeneTestUtils.createTransExons;
import static com.hartwig.hmftools.common.test.GeneTestUtils.generateExonStarts;
import static com.hartwig.hmftools.common.test.GeneTestUtils.generateTransName;
import static com.hartwig.hmftools.common.gene.TranscriptProteinData.BIOTYPE_PROTEIN_CODING;
import static com.hartwig.hmftools.common.fusion.FusionCommon.NEG_STRAND;
import static com.hartwig.hmftools.common.fusion.FusionCommon.POS_STRAND;
import static com.hartwig.hmftools.common.fusion.KnownFusionType.EXON_DEL_DUP;
import static com.hartwig.hmftools.common.fusion.KnownFusionType.IG_KNOWN_PAIR;
import static com.hartwig.hmftools.common.fusion.KnownFusionType.IG_PROMISCUOUS;
import static com.hartwig.hmftools.common.fusion.KnownFusionType.KNOWN_PAIR;
import static com.hartwig.hmftools.common.fusion.KnownFusionType.PROMISCUOUS_3;
import static com.hartwig.hmftools.common.genome.region.Orientation.ORIENT_REV;
import static com.hartwig.hmftools.common.genome.region.Orientation.ORIENT_FWD;
import static com.hartwig.hmftools.linx.fusion.FusionConstants.PRE_GENE_PROMOTOR_DISTANCE;
import static com.hartwig.hmftools.linx.gene.BreakendGenePrep.findGeneAnnotationsBySv;
import static com.hartwig.hmftools.linx.utils.GeneTestUtils.CHR_1;
import static com.hartwig.hmftools.linx.utils.GeneTestUtils.CHR_2;
import static com.hartwig.hmftools.linx.utils.GeneTestUtils.GENE_NAME_1;
import static com.hartwig.hmftools.linx.utils.GeneTestUtils.GENE_NAME_2;
import static com.hartwig.hmftools.linx.utils.GeneTestUtils.GENE_NAME_3;
import static com.hartwig.hmftools.linx.utils.GeneTestUtils.TRANS_1;
import static com.hartwig.hmftools.linx.utils.GeneTestUtils.addTestGenes;
import static com.hartwig.hmftools.linx.utils.GeneTestUtils.addTestTranscripts;
import static com.hartwig.hmftools.linx.utils.SampleDataLoader.setSvGeneData;
import static com.hartwig.hmftools.linx.utils.SvTestUtils.createBnd;
import static com.hartwig.hmftools.linx.utils.SvTestUtils.createDel;
import static com.hartwig.hmftools.linx.utils.SvTestUtils.createInf;
import static com.hartwig.hmftools.linx.utils.SvTestUtils.createSgl;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertNotNull;
import static org.junit.Assert.assertTrue;

import java.util.List;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.ensemblcache.EnsemblDataCache;
import com.hartwig.hmftools.common.gene.GeneData;
import com.hartwig.hmftools.common.gene.TranscriptData;
import com.hartwig.hmftools.linx.gene.BreakendGeneData;
import com.hartwig.hmftools.common.fusion.KnownFusionData;
import com.hartwig.hmftools.linx.types.SglMapping;
import com.hartwig.hmftools.linx.types.SvCluster;
import com.hartwig.hmftools.linx.types.SvVarData;
import com.hartwig.hmftools.linx.utils.LinxTester;

import org.junit.Test;

public class SpecialFusionsTest
{
    @Test
    public void testSameGeneFusions()
    {
        LinxTester tester = new LinxTester();

        EnsemblDataCache geneTransCache = createGeneDataCache();

        tester.initialiseFusions(geneTransCache);

        addTestGenes(geneTransCache);
        addTestTranscripts(geneTransCache);

        // mark as 3' promiscuous
        tester.FusionAnalyser.getFusionFinder().getKnownFusionCache()
                .addData(new KnownFusionData(PROMISCUOUS_3, "", GENE_NAME_1, "", ""));

        List<BreakendGeneData> upGenes = Lists.newArrayList();
        int upPos = 1950;
        upGenes.addAll(findGeneAnnotationsBySv(geneTransCache, 0, false, CHR_1, upPos, ORIENT_FWD, PRE_GENE_PROMOTOR_DISTANCE));
        upGenes.get(0).setPositionalData(CHR_1, upPos, ORIENT_FWD);

        // add downstream breakends
        List<BreakendGeneData> downGenes = Lists.newArrayList();
        int downPos = 1350;
        downGenes.addAll(findGeneAnnotationsBySv(geneTransCache, 0, true, CHR_1, downPos, ORIENT_REV, PRE_GENE_PROMOTOR_DISTANCE));
        downGenes.get(0).setPositionalData(CHR_1, downPos, ORIENT_REV);

        tester.FusionAnalyser.getFusionFinder().setFusionParams(true, true, false);

        List<GeneFusion> fusions = tester.FusionAnalyser.getFusionFinder().findFusions(upGenes, downGenes);

        assertEquals(1, fusions.size());
        final GeneFusion fusion = fusions.get(0);

        assertEquals(upPos, fusion.upstreamTrans().breakendGeneData().position());
        assertEquals(downPos, fusion.downstreamTrans().breakendGeneData().position());
        assertEquals(0, fusion.getExonsSkipped(false));
        assertEquals(0, fusion.getExonsSkipped(false));
        assertTrue(!fusion.reportable());
    }

    @Test
    public void testExonDelDupFusion()
    {
        LinxTester tester = new LinxTester();

        EnsemblDataCache geneTransCache = createGeneDataCache();

        tester.initialiseFusions(geneTransCache);

        addTestGenes(geneTransCache);
        addTestTranscripts(geneTransCache);

        String transName = generateTransName(TRANS_1);

        // String knownDelRegion = String.format("%s;%d;%d;%d;%d", transName, 2, 3, 5, 6);

        KnownFusionData kfData = new KnownFusionData(EXON_DEL_DUP, GENE_NAME_1, GENE_NAME_1, "", "");
        kfData.setKnownExonData(transName, "2;3", "5;6");
        tester.FusionAnalyser.getFusionFinder().getKnownFusionCache().addData(kfData);

        tester.FusionAnalyser.getFusionFinder().setFusionParams(true, true, false);

        // first DEL doesn't delete a known region even though it's phased
        List<BreakendGeneData> upGenes = Lists.newArrayList();
        int upPos = 1550;
        upGenes.addAll(findGeneAnnotationsBySv(geneTransCache, 0, true, CHR_1, upPos, ORIENT_FWD, PRE_GENE_PROMOTOR_DISTANCE));
        upGenes.get(0).setPositionalData(CHR_1, upPos, ORIENT_FWD);

        // add downstream breakends
        List<BreakendGeneData> downGenes = Lists.newArrayList();

        int downPos = 2150;
        downGenes.addAll(findGeneAnnotationsBySv(geneTransCache, 0, false, CHR_1, downPos, ORIENT_REV, PRE_GENE_PROMOTOR_DISTANCE));
        downGenes.get(0).setPositionalData(CHR_1, downPos, ORIENT_REV);

        List<GeneFusion> fusions = tester.FusionAnalyser.getFusionFinder().findFusions(upGenes, downGenes);

        // second one does
        upGenes.clear();
        downGenes.clear();
        upPos = 1350;
        upGenes.addAll(findGeneAnnotationsBySv(geneTransCache, 1, true, CHR_1, upPos, ORIENT_FWD, PRE_GENE_PROMOTOR_DISTANCE));
        upGenes.get(0).setPositionalData(CHR_1, upPos, ORIENT_FWD);

        downPos = 1950;
        downGenes.addAll(findGeneAnnotationsBySv(geneTransCache, 1, false, CHR_1, downPos, ORIENT_REV, PRE_GENE_PROMOTOR_DISTANCE));
        downGenes.get(0).setPositionalData(CHR_1, downPos, ORIENT_REV);

        fusions.addAll(tester.FusionAnalyser.getFusionFinder().findFusions(upGenes, downGenes));

        assertEquals(1, fusions.size());
        GeneFusion fusion = tester.FusionAnalyser.getFusionFinder().findTopReportableFusion(fusions);
        assertTrue(fusion != null);
        assertEquals(EXON_DEL_DUP, fusion.knownType());

        // the selected fusion is the longest for coding bases and without any exon skipping
        assertEquals(upPos, fusion.upstreamTrans().breakendGeneData().position());
        assertEquals(downPos, fusion.downstreamTrans().breakendGeneData().position());
        assertEquals(0, fusion.getExonsSkipped(true));
        assertEquals(0, fusion.getExonsSkipped(false));
        assertTrue(fusion.reportable());

        // test again with exon skipping as long as the skipping occurs within the bounds specified for the exon DEL-DUP

        List<GeneData> geneList = Lists.newArrayList();
        final String geneId = "ENSG0010";
        int geneStart = 30000;
        geneList.add(createEnsemblGeneData(geneId, geneId, CHR_1, POS_STRAND, geneStart, geneStart + 15 * 200));
        addGeneData(geneTransCache, CHR_1, geneList);

        int transId = 10;
        TranscriptData transData = createTransExons(
                geneId, transId, POS_STRAND, generateExonStarts(30000, 15, 100, 100),
                100, geneStart + 250, geneStart + 2850, true, BIOTYPE_PROTEIN_CODING);

        addTransExonData(geneTransCache, geneId, Lists.newArrayList(transData));

        transName = generateTransName(transId);

        kfData = new KnownFusionData(EXON_DEL_DUP, geneId, geneId, "", "");
        kfData.setKnownExonData(transName, "2;5", "10;14");

        tester.FusionAnalyser.getFusionFinder().getKnownFusionCache().addData(kfData);

        upGenes.clear();
        downGenes.clear();
        upPos = geneStart + 350;
        upGenes.addAll(findGeneAnnotationsBySv(geneTransCache, 1, true, CHR_1, upPos, ORIENT_FWD, PRE_GENE_PROMOTOR_DISTANCE));
        upGenes.get(0).setPositionalData(CHR_1, upPos, ORIENT_FWD);

        downPos = geneStart + 1950;
        downGenes.addAll(findGeneAnnotationsBySv(geneTransCache, 1, false, CHR_1, downPos, ORIENT_REV, PRE_GENE_PROMOTOR_DISTANCE));
        downGenes.get(0).setPositionalData(CHR_1, downPos, ORIENT_REV);

        fusions.clear();
        fusions.addAll(tester.FusionAnalyser.getFusionFinder().findFusions(upGenes, downGenes));

        assertEquals(1, fusions.size());
        fusion = tester.FusionAnalyser.getFusionFinder().findTopReportableFusion(fusions);
        assertTrue(fusion != null);
        assertEquals(EXON_DEL_DUP, fusion.knownType());

        // the selected fusion is the longest for coding bases and without any exon skipping
        assertEquals(upPos, fusion.upstreamTrans().breakendGeneData().position());
        assertEquals(downPos, fusion.downstreamTrans().breakendGeneData().position());
        assertEquals(0, fusion.getExonsSkipped(true));
        assertEquals(1, fusion.getExonsSkipped(false));
        assertEquals(2, fusion.getFusedExon(true));
        assertEquals(12, fusion.getFusedExon(false));
        assertTrue(fusion.reportable());
    }

    @Test
    public void testIgRegionFusion()
    {
        LinxTester tester = new LinxTester();

        EnsemblDataCache geneTransCache = createGeneDataCache();

        tester.initialiseFusions(geneTransCache);

        String geneName = "IGH";
        String geneId1 = "ENSG0001";
        String chromosome = "1";

        String geneName2 = "GENE2";
        String geneId2 = "ENSG0002";

        String geneName3 = "GENE3";
        String geneId3 = "ENSG0003";

        byte strand = POS_STRAND;

        List<GeneData> geneList = Lists.newArrayList();
        geneList.add(createEnsemblGeneData(geneId1, geneName, chromosome, strand, 100, 1500));
        geneList.add(createEnsemblGeneData(geneId2, geneName2, chromosome, strand, 10000, 11500));
        geneList.add(createEnsemblGeneData(geneId3, geneName3, chromosome, strand, 20000, 21500));

        addGeneData(geneTransCache, chromosome, geneList);

        List<TranscriptData> transDataList = Lists.newArrayList();

        int transId = 1;

        int[] exonStarts = generateExonStarts(100, 7, 98, 100);

        TranscriptData transData = createTransExons(geneId1, transId++, strand, exonStarts, 98, 300, 1398, true, "");
        transDataList.add(transData);

        addTransExonData(geneTransCache, geneId1, transDataList);

        exonStarts = generateExonStarts(10000, exonStarts.length, 98, 100);
        transData = createTransExons(geneId2, transId++, strand, exonStarts, 98, 10300, 11390, true, "");
        transDataList = Lists.newArrayList(transData);

        addTransExonData(geneTransCache, geneId2, transDataList);

        exonStarts = generateExonStarts(20000, exonStarts.length, 98, 100);
        transData = createTransExons(geneId3, transId++, strand, exonStarts, 98, 20300, 21390, true, "");
        transDataList = Lists.newArrayList(transData);

        addTransExonData(geneTransCache, geneId3, transDataList);

        final String igRegion = String.format("IG_RANGE=%d;%s;%d;%d", NEG_STRAND, chromosome, 50, 2000);

        KnownFusionData kfData = new KnownFusionData(IG_KNOWN_PAIR, geneName, geneName2, "", "");
        kfData.applyOverrides(igRegion);

        tester.FusionAnalyser.getFusionFinder().getKnownFusionCache().addData(kfData);

        kfData = new KnownFusionData(IG_PROMISCUOUS, geneName, "", "", "");
        kfData.applyOverrides(igRegion);

        tester.FusionAnalyser.getFusionFinder().getKnownFusionCache().addData(kfData);

        tester.FusionAnalyser.getFusionFinder().setFusionParams(true, true, false);

        // a DEL linking the 2 regions
        List<BreakendGeneData> upGenes = Lists.newArrayList();
        upGenes.addAll(findGeneAnnotationsBySv(geneTransCache, 0, true, chromosome, 200, ORIENT_REV, 1000));
        upGenes.get(0).setPositionalData(chromosome, 200, ORIENT_REV);

        // add downstream breakends
        List<BreakendGeneData> downGenes = Lists.newArrayList();

        downGenes.addAll(findGeneAnnotationsBySv(geneTransCache, 0, false, chromosome, 9500, ORIENT_REV, 1000));
        downGenes.get(0).setPositionalData(chromosome, 9500, ORIENT_REV);

        List<GeneFusion> fusions = tester.FusionAnalyser.getFusionFinder().findFusions(upGenes, downGenes);

        upGenes.clear();
        upGenes.addAll(findGeneAnnotationsBySv(geneTransCache, 1, true, chromosome, 200, ORIENT_REV, 1000));
        upGenes.get(0).setPositionalData(chromosome, 500, ORIENT_REV);

        downGenes.clear();
        downGenes.addAll(findGeneAnnotationsBySv(geneTransCache, 1, false, chromosome, 19500, ORIENT_REV, 1000));
        downGenes.get(0).setPositionalData(chromosome, 20100, ORIENT_REV);

        fusions.addAll(tester.FusionAnalyser.getFusionFinder().findFusions(upGenes, downGenes));

        assertEquals(2, fusions.size());

        GeneFusion fusion = tester.FusionAnalyser.getFusionFinder().findTopReportableFusion(fusions);
        assertTrue(fusion != null);
        assertEquals(IG_KNOWN_PAIR, fusion.knownType());

        // the selected fusion is the longest for coding bases and without any exon skipping
        assertEquals(200, fusion.upstreamTrans().breakendGeneData().position());
        assertEquals(9500, fusion.downstreamTrans().breakendGeneData().position());
        assertTrue(fusion.reportable());

        fusion = fusions.stream().filter(x -> x.knownType() == IG_PROMISCUOUS).findFirst().orElse(null);
        assertTrue(fusion != null);

        // the selected fusion is the longest for coding bases and without any exon skipping
        assertEquals(500, fusion.upstreamTrans().breakendGeneData().position());
        assertEquals(20100, fusion.downstreamTrans().breakendGeneData().position());
        assertTrue(!fusion.reportable());
    }

    @Test
    public void testSingleBreakendFusions()
    {
        LinxTester tester = new LinxTester();

        EnsemblDataCache geneTransCache = createGeneDataCache();

        tester.initialiseFusions(geneTransCache);

        addTestGenes(geneTransCache);
        addTestTranscripts(geneTransCache);

        PRE_GENE_PROMOTOR_DISTANCE = 200;

        // set known fusion gene for the SGL breakend

        KnownFusionData kfData = new KnownFusionData(KNOWN_PAIR, GENE_NAME_1, GENE_NAME_2, "", "");
        tester.FusionAnalyser.getFusionFinder().getKnownFusionCache().addData(kfData);

        // test 1: a SGL by itself
        int varId = 1;

        SvVarData sgl1 = createSgl(varId++, CHR_1, 1150, ORIENT_FWD);

        final String altMapping = CHR_1 + ":" + String.valueOf(10150) + "|" + "+" + "|" + "10M" + "|19";
        sgl1.getSglMappings().add(SglMapping.from(altMapping, ORIENT_FWD));

        tester.AllVariants.add(sgl1);

        tester.preClusteringInit();
        tester.Analyser.clusterAndAnalyse();

        setSvGeneData(tester.AllVariants, geneTransCache, false);
        tester.FusionAnalyser.annotateTranscripts(tester.AllVariants, false);

        tester.FusionAnalyser.run(tester.SampleId, tester.AllVariants, tester.getClusters(), tester.Analyser.getState().getChrBreakendMap());

        assertEquals(1, tester.FusionAnalyser.getFusions().size());

        GeneFusion fusion = tester.FusionAnalyser.getFusions().get(0);
        assertEquals(sgl1.id(), fusion.upstreamTrans().breakendGeneData().varId());
        assertEquals(sgl1.id(), fusion.downstreamTrans().breakendGeneData().varId());
        assertEquals(10150, fusion.downstreamTrans().breakendGeneData().position());

        tester.clearClustersAndSVs();

        // test 2: test a chain with a SGL at the start
        sgl1 = createSgl(varId++, CHR_1, 10150, ORIENT_REV);

        sgl1.getSglMappings().add(new SglMapping(CHR_1, 1150, ORIENT_FWD, "", 1));

        SvVarData var1 = createDel(varId++, CHR_1, 14000,15000);

        tester.AllVariants.add(sgl1);
        tester.AllVariants.add(var1);

        tester.preClusteringInit();
        tester.Analyser.clusterAndAnalyse();

        assertEquals(1, tester.Analyser.getClusters().size());
        SvCluster cluster = tester.Analyser.getClusters().get(0);

        assertEquals(1, cluster.getChains().size());

        setSvGeneData(tester.AllVariants, geneTransCache, true);
        tester.FusionAnalyser.annotateTranscripts(tester.AllVariants, true);

        tester.FusionAnalyser.run(tester.SampleId, tester.AllVariants, tester.getClusters(), tester.Analyser.getState().getChrBreakendMap());

        assertEquals(1, tester.FusionAnalyser.getFusions().size());

        fusion = tester.FusionAnalyser.getFusions().get(0);
        assertEquals(sgl1.id(), fusion.upstreamTrans().breakendGeneData().varId());
        assertEquals(sgl1.id(), fusion.downstreamTrans().breakendGeneData().varId());

        // test 3: a chain with the SGL breaking going outside the gene into a shard
        tester.clearClustersAndSVs();

        sgl1 = createSgl(varId++, CHR_2, 100, ORIENT_REV);

        sgl1.getSglMappings().add(new SglMapping(CHR_1, 1150, ORIENT_FWD, "", 1));

        var1 = createBnd(varId++, CHR_1, 10150, ORIENT_REV, CHR_2, 200, ORIENT_FWD);
        SvVarData var2 = createBnd(varId++, CHR_1, 500, ORIENT_FWD, CHR_2, 500, ORIENT_REV); // to assist clustering only

        tester.AllVariants.add(sgl1);
        tester.AllVariants.add(var1);
        tester.AllVariants.add(var2);

        tester.preClusteringInit();
        tester.Analyser.clusterAndAnalyse();

        assertEquals(1, tester.Analyser.getClusters().size());
        cluster = tester.Analyser.getClusters().get(0);

        assertEquals(1, cluster.getChains().size());

        setSvGeneData(tester.AllVariants, geneTransCache, true);
        tester.FusionAnalyser.annotateTranscripts(tester.AllVariants, true);

        tester.FusionAnalyser.run(tester.SampleId, tester.AllVariants, tester.getClusters(), tester.Analyser.getState().getChrBreakendMap());

        assertEquals(1, tester.FusionAnalyser.getFusions().size());

        fusion = tester.FusionAnalyser.getFusions().get(0);
        assertEquals(sgl1.id(), fusion.upstreamTrans().breakendGeneData().varId());
        assertEquals(var1.id(), fusion.downstreamTrans().breakendGeneData().varId());

        // test 4: 2 SGLs at the start and end of a fusion
        tester.clearClustersAndSVs();

        sgl1 = createSgl(varId++, CHR_2, 100, ORIENT_REV);
        sgl1.getSglMappings().add(new SglMapping(CHR_1, 1150, ORIENT_FWD, "", 1));

        var1 = createBnd(varId++, CHR_2, 200, ORIENT_FWD, "3", 200, ORIENT_FWD);

        SvVarData sgl2 = createSgl(varId++, "3", 100, ORIENT_REV);
        sgl2.getSglMappings().add(new SglMapping(CHR_1, 10150, ORIENT_REV, "", 1));

        tester.AllVariants.add(sgl1);
        tester.AllVariants.add(sgl2);
        tester.AllVariants.add(var1);

        tester.preClusteringInit();
        tester.Analyser.clusterAndAnalyse();

        assertEquals(1, tester.Analyser.getClusters().size());
        cluster = tester.Analyser.getClusters().get(0);

        assertEquals(1, cluster.getChains().size());
        assertEquals(3, cluster.getChains().get(0).getSvCount());

        setSvGeneData(tester.AllVariants, geneTransCache, true);
        tester.FusionAnalyser.annotateTranscripts(tester.AllVariants, true);

        tester.FusionAnalyser.run(tester.SampleId, tester.AllVariants, tester.getClusters(), tester.Analyser.getState().getChrBreakendMap());

        assertEquals(1, tester.FusionAnalyser.getFusions().size());

        fusion = tester.FusionAnalyser.getFusions().get(0);
        assertEquals(sgl1.id(), fusion.upstreamTrans().breakendGeneData().varId());
        assertEquals(sgl2.id(), fusion.downstreamTrans().breakendGeneData().varId());
    }

    @Test
    public void testSglMappedInfFusions()
    {
        LinxTester tester = new LinxTester();

        EnsemblDataCache geneTransCache = createGeneDataCache();

        tester.initialiseFusions(geneTransCache);

        addTestGenes(geneTransCache);
        addTestTranscripts(geneTransCache);

        KnownFusionData kfData = new KnownFusionData(KNOWN_PAIR, GENE_NAME_1, GENE_NAME_2, "", "");
        tester.FusionAnalyser.getFusionFinder().getKnownFusionCache().addData(kfData);

        int varId = 1;

        SvVarData sgl = createSgl(varId++, CHR_1, 1150, ORIENT_FWD);
        SvVarData inf = createInf(varId++, CHR_2, 100050, ORIENT_REV);

        final String altMapping = CHR_2 + ":" + String.valueOf(100000) + "|" + "+" + "|" + "10M" + "|19";
        sgl.getSglMappings().add(SglMapping.from(altMapping, ORIENT_FWD));

        SvVarData bnd = createBnd(varId++, CHR_1, 10020, ORIENT_REV, CHR_2, 100200, ORIENT_FWD);

        tester.AllVariants.add(sgl);
        tester.AllVariants.add(inf);
        tester.AllVariants.add(bnd);

        tester.preClusteringInit();

        assertEquals(4, tester.AllVariants.size());
        SvVarData sglInf = tester.AllVariants.get(3);

        tester.Analyser.clusterAndAnalyse();

        setSvGeneData(tester.AllVariants, geneTransCache, false);
        tester.FusionAnalyser.annotateTranscripts(tester.AllVariants, false);

        tester.FusionAnalyser.run(tester.SampleId, tester.AllVariants, tester.getClusters(), tester.Analyser.getState().getChrBreakendMap());

        assertEquals(1, tester.FusionAnalyser.getFusions().size());

        GeneFusion fusion = tester.FusionAnalyser.getFusions().get(0);
        assertEquals(sglInf.id(), fusion.upstreamTrans().breakendGeneData().varId());
        assertEquals(bnd.id(), fusion.downstreamTrans().breakendGeneData().varId());
    }

    @Test
    public void testKnownUnmappableFusions()
    {
        LinxTester tester = new LinxTester();

        EnsemblDataCache geneTransCache = createGeneDataCache();

        tester.initialiseFusions(geneTransCache);

        addTestGenes(geneTransCache);
        addTestTranscripts(geneTransCache);

        String igGeneName = "IGH";
        String igGeneId = "ENSG0IGH";
        String CHR_4 = "4";

        byte strand = POS_STRAND;

        List<GeneData> geneList = Lists.newArrayList();
        geneList.add(createEnsemblGeneData(igGeneId, igGeneName, CHR_4, strand, 100, 1500));

        addGeneData(geneTransCache, CHR_4, geneList);

        List<TranscriptData> transDataList = Lists.newArrayList();

        int transId = 1;

        int[] exonStarts = generateExonStarts(100, 7, 100, 100);

        int codingStart = 100;
        int codingEnd = 799;
        TranscriptData transData = createTransExons(igGeneId, transId++, strand, exonStarts, 100, codingStart, codingEnd, true, "");
        transDataList.add(transData);

        addTransExonData(geneTransCache, igGeneId, transDataList);

        // first fusion goes from

        final String igRegionStr = String.format("IG_RANGE=%d;%s;%d;%d", POS_STRAND, CHR_4, 50, 2000);

        // set known fusion gene for the SGL breakend
        // the 3' gene's alt mapping is about expanding its range beyond the Ensembl-defined one
        final String threeAltMapping = String.format("ALTS=ALT;1;20000;25000;ALT;GS;50000;60000");

        KnownFusionData kfData = new KnownFusionData(KNOWN_PAIR, GENE_NAME_1, GENE_NAME_3, "", "");
        kfData.applyOverrides(threeAltMapping);
        tester.FusionAnalyser.getFusionFinder().getKnownFusionCache().addData(kfData);

        kfData = new KnownFusionData(IG_KNOWN_PAIR, igGeneName, GENE_NAME_3, "", "");
        kfData.applyOverrides(String.format("%s %s", igRegionStr, threeAltMapping));

        tester.FusionAnalyser.getFusionFinder().getKnownFusionCache().addData(kfData);

        tester.FusionAnalyser.getSpecialFusions().cacheSpecialFusionGenes();

        // first test a fusion which uses the expanded 3' alt mapping to gene 3
        int varId = 1;
        SvVarData sgl1 = createSgl(varId++, CHR_1, 1150, ORIENT_FWD);
        sgl1.getSglMappings().add(new SglMapping(CHR_1, 24000, ORIENT_FWD, "", 1));

        tester.AllVariants.add(sgl1);

        tester.preClusteringInit();
        tester.Analyser.clusterAndAnalyse();

        setSvGeneData(tester.AllVariants, geneTransCache, false);
        tester.FusionAnalyser.annotateTranscripts(tester.AllVariants, false);

        tester.FusionAnalyser.run(tester.SampleId, tester.AllVariants, tester.getClusters(), tester.Analyser.getState().getChrBreakendMap());

        assertEquals(1, tester.FusionAnalyser.getFusions().size());

        GeneFusion fusion = tester.FusionAnalyser.getFusions().get(0);
        assertNotNull(fusion);
        assertEquals(sgl1.id(), fusion.upstreamTrans().breakendGeneData().varId());
        assertEquals(sgl1.id(), fusion.downstreamTrans().breakendGeneData().varId());
        assertTrue(fusion.reportable());

        assertEquals(KNOWN_PAIR, fusion.knownType());

        // try again but to an alternate mapping location
        tester.clearClustersAndSVs();
        sgl1 = createSgl(varId++, CHR_1, 1150, ORIENT_FWD);
        sgl1.getSglMappings().add(new SglMapping("GS", 55000, ORIENT_FWD, "", 1));

        tester.AllVariants.add(sgl1);

        tester.preClusteringInit();
        tester.Analyser.clusterAndAnalyse();

        setSvGeneData(tester.AllVariants, geneTransCache, false);
        tester.FusionAnalyser.annotateTranscripts(tester.AllVariants, false);

        tester.FusionAnalyser.run(tester.SampleId, tester.AllVariants, tester.getClusters(), tester.Analyser.getState().getChrBreakendMap());

        assertEquals(1, tester.FusionAnalyser.getFusions().size());

        fusion = tester.FusionAnalyser.getFusions().get(0);
        assertEquals(sgl1.id(), fusion.upstreamTrans().breakendGeneData().varId());
        assertEquals(sgl1.id(), fusion.downstreamTrans().breakendGeneData().varId());
        assertTrue(fusion.reportable());

        assertEquals(KNOWN_PAIR, fusion.knownType());

        // test again but using an IG region for the 5' gene
        tester.clearClustersAndSVs();
        sgl1 = createSgl(varId++, CHR_4, 1000, ORIENT_FWD);
        sgl1.getSglMappings().add(new SglMapping("GS", 55000, ORIENT_FWD, "", 1));

        tester.AllVariants.add(sgl1);

        tester.preClusteringInit();
        tester.Analyser.clusterAndAnalyse();

        setSvGeneData(tester.AllVariants, geneTransCache, false);
        tester.FusionAnalyser.annotateTranscripts(tester.AllVariants, false);

        tester.FusionAnalyser.run(tester.SampleId, tester.AllVariants, tester.getClusters(), tester.Analyser.getState().getChrBreakendMap());

        assertEquals(1, tester.FusionAnalyser.getFusions().size());

        fusion = tester.FusionAnalyser.getFusions().get(0);
        assertTrue(fusion.reportable());
        assertEquals(IG_KNOWN_PAIR, fusion.knownType());
        assertEquals(igGeneName, fusion.upstreamTrans().geneName());
        assertEquals(GENE_NAME_3, fusion.downstreamTrans().geneName());
    }
}
