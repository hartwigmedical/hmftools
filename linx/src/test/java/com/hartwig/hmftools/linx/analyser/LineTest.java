package com.hartwig.hmftools.linx.analyser;

import static com.hartwig.hmftools.common.utils.sv.SvCommonUtils.NEG_ORIENT;
import static com.hartwig.hmftools.common.utils.sv.SvCommonUtils.POS_ORIENT;
import static com.hartwig.hmftools.common.sv.StructuralVariantType.BND;
import static com.hartwig.hmftools.common.sv.StructuralVariantType.DEL;
import static com.hartwig.hmftools.common.sv.StructuralVariantType.INV;
import static com.hartwig.hmftools.common.sv.StructuralVariantType.SGL;
import static com.hartwig.hmftools.linx.annotators.LineClusterState.containsMotif;
import static com.hartwig.hmftools.linx.annotators.LineClusterState.hasLinePolyAorTMotif;
import static com.hartwig.hmftools.linx.annotators.LineElementAnnotator.POLY_A_MOTIF;
import static com.hartwig.hmftools.linx.annotators.LineElementAnnotator.POLY_T_MOTIF;
import static com.hartwig.hmftools.linx.annotators.LineElementType.KNOWN;
import static com.hartwig.hmftools.linx.annotators.LineElementType.SUSPECT;
import static com.hartwig.hmftools.linx.utils.GeneTestUtils.CHR_2;
import static com.hartwig.hmftools.linx.utils.SvTestUtils.createBnd;
import static com.hartwig.hmftools.linx.utils.SvTestUtils.createInf;
import static com.hartwig.hmftools.linx.utils.SvTestUtils.createInv;
import static com.hartwig.hmftools.linx.utils.SvTestUtils.createSgl;
import static com.hartwig.hmftools.linx.utils.SvTestUtils.createSv;
import static com.hartwig.hmftools.linx.utils.SvTestUtils.createTestSv;

import static org.junit.Assert.assertFalse;

import static junit.framework.TestCase.assertEquals;
import static junit.framework.TestCase.assertTrue;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.linx.annotators.LineElementAnnotator;
import com.hartwig.hmftools.linx.types.SglMapping;
import com.hartwig.hmftools.linx.types.SvCluster;
import com.hartwig.hmftools.linx.types.SvVarData;
import com.hartwig.hmftools.linx.utils.LinxTester;

import org.junit.Test;

public class LineTest
{
    @Test
    public void testKnownLineMarking()
    {
        LinxTester tester = new LinxTester();

        // scenario 1:
        // 2 BNDs within 5KB and 1 have poly A/T and not forming a DB
        // if other breakends form a DB, don't mark
        SvVarData bnd1 = createBnd(tester.nextVarId(), "1", 100, -1, "2", 100, 1);
        SvVarData bnd2 = createBnd(tester.nextVarId(), "1", 200, 1, "2", 200, -1);

        // don't need to be marked as LINE
        SvVarData sgl = createSgl(tester.nextVarId(), "1", 500, 1);
        SvVarData inf = createInf(tester.nextVarId(), "1", 1500, 1);

        bnd1.addLineElement(KNOWN, true);
        bnd2.addLineElement(KNOWN, true);

        SvCluster cluster = new SvCluster(0);
        cluster.addVariant(bnd1);
        cluster.addVariant(bnd2);
        cluster.addVariant(sgl);
        cluster.addVariant(inf);

        tester.addClusterAndSVs(cluster);

        int proximity = 5000;
        LineElementAnnotator leAnnotator = new LineElementAnnotator(proximity);

        leAnnotator.markLineCluster(cluster);

        assertTrue(cluster.hasLinkingLineElements());
    }
    @Test
    public void testLineMotifs()
    {
        String shortBases = "GCTG";
        String randomBases = "GCTGTGAGCTGATCG";

        String insSequence = POLY_A_MOTIF.substring(0, 8) + shortBases + POLY_A_MOTIF.substring(0, 8);
        assertTrue(containsMotif(insSequence, POLY_A_MOTIF));

        insSequence = POLY_A_MOTIF.substring(0, 7) + shortBases + POLY_A_MOTIF.substring(0, 8);
        assertFalse(containsMotif(insSequence, POLY_A_MOTIF));

        insSequence = shortBases + POLY_A_MOTIF + randomBases;

        assertFalse(hasLinePolyAorTMotif(shortBases, POS_ORIENT, true, false));
        assertFalse(hasLinePolyAorTMotif(randomBases, POS_ORIENT, true, false));

        assertTrue(hasLinePolyAorTMotif(insSequence, POS_ORIENT, true, false));
        assertFalse(hasLinePolyAorTMotif(insSequence, POS_ORIENT, false, false));
        assertFalse(hasLinePolyAorTMotif(insSequence, NEG_ORIENT, true, false));
        assertFalse(hasLinePolyAorTMotif(insSequence, NEG_ORIENT, false, false));

        insSequence = shortBases + POLY_T_MOTIF + randomBases;

        assertFalse(hasLinePolyAorTMotif(insSequence, POS_ORIENT, true, false));
        assertTrue(hasLinePolyAorTMotif(insSequence, POS_ORIENT, false, false));
        assertFalse(hasLinePolyAorTMotif(insSequence, NEG_ORIENT, true, false));
        assertFalse(hasLinePolyAorTMotif(insSequence, NEG_ORIENT, false, false));

        insSequence = randomBases + POLY_T_MOTIF + shortBases;

        assertFalse(hasLinePolyAorTMotif(insSequence, POS_ORIENT, true, false));
        assertFalse(hasLinePolyAorTMotif(insSequence, POS_ORIENT, false, false));
        assertTrue(hasLinePolyAorTMotif(insSequence, NEG_ORIENT, true, false));
        assertFalse(hasLinePolyAorTMotif(insSequence, NEG_ORIENT, false, false));

        insSequence = randomBases + POLY_A_MOTIF + shortBases;

        assertFalse(hasLinePolyAorTMotif(insSequence, POS_ORIENT, true, false));
        assertFalse(hasLinePolyAorTMotif(insSequence, POS_ORIENT, false, false));
        assertFalse(hasLinePolyAorTMotif(insSequence, NEG_ORIENT, true, false));
        assertTrue(hasLinePolyAorTMotif(insSequence, NEG_ORIENT, false, false));

        // test allowing the sequence at either end
        assertFalse(hasLinePolyAorTMotif(insSequence, POS_ORIENT, true, false));
        assertTrue(hasLinePolyAorTMotif(insSequence, POS_ORIENT, true, true));

        insSequence = shortBases + POLY_T_MOTIF + randomBases;
        assertFalse(hasLinePolyAorTMotif(insSequence, NEG_ORIENT, true, false));
        assertTrue(hasLinePolyAorTMotif(insSequence, NEG_ORIENT, true, true));
    }

    @Test
    public void testSuspectBndLineMarking()
    {
        LinxTester tester = new LinxTester();

        String shortBases = "GCTG";
        String randomBases = "GCTGTGAGCTGATCG";
        String polyAMotif = shortBases + POLY_A_MOTIF + randomBases;
        String polyTMotif = randomBases + POLY_T_MOTIF + shortBases;
        String insertPolyAMotif = randomBases + POLY_A_MOTIF + shortBases;

        // scenario 1:
        // 2 BNDs within 5KB and 1 have poly A/T and not forming a DB with other breakends forming a short DB
        SvVarData bnd1 = createSv(tester.nextVarId(), "1", "2", 100, 150,  -1, 1, BND, polyTMotif);
        SvVarData bnd2 = createSv(tester.nextVarId(), "1", "2", 200, 170,  1, -1, BND, "");

        tester.addAndCluster(bnd1, bnd2);

        // in this case because the DB isn't short on the originating arm, both ends are marked as suspect
        assertTrue(bnd1.hasLineElement(SUSPECT, true));
        assertTrue(bnd2.hasLineElement(SUSPECT, true));
        assertFalse(bnd1.isLineElement(false));
        assertFalse(bnd2.isLineElement(false));

        SvCluster cluster = tester.Analyser.getClusters().get(0);
        assertTrue(cluster.hasLinkingLineElements());

        // 2 BNDs both with poly A/T
        bnd1 = createSv(tester.nextVarId(), "1", "2", 100, 100,  -1, 1, BND, polyTMotif);
        bnd2 = createSv(tester.nextVarId(), "1", "3", 200, 110,  1, -1, BND, polyAMotif);

        tester.addAndCluster(bnd1, bnd2);

        assertTrue(bnd1.hasLineElement(SUSPECT, true));
        assertTrue(bnd2.hasLineElement(SUSPECT, true));
        assertFalse(bnd1.isLineElement(false));
        assertFalse(bnd2.isLineElement(false));

        cluster = tester.Analyser.getClusters().get(0);
        assertTrue(cluster.hasLinkingLineElements());

        // now a BND with a remote SGL in a short DB
        bnd1 = createSv(tester.nextVarId(), "1", "2", 100, 100,  -1, 1, BND, polyTMotif);
        SvVarData sgl = createSv(tester.nextVarId(), "2", "", 90, -1,  -1, 0, SGL, insertPolyAMotif);

        tester.addAndCluster(bnd1, sgl);

        assertTrue(bnd1.hasLineElement(SUSPECT, true));
        assertFalse(sgl.isLineElement(true));
        assertFalse(bnd1.isLineElement(false));

        cluster = tester.Analyser.getClusters().get(0);
        assertTrue(cluster.hasLinkingLineElements());

        // now with BNDs going to different arms
        bnd1 = createSv(tester.nextVarId(), "1", "2", 100, 100,  -1, 1, BND, polyTMotif);
        bnd2 = createSv(tester.nextVarId(), "1", "3", 200, 110,  1, -1, BND, "");

        tester.addAndCluster(bnd1, bnd2);

        assertTrue(bnd1.hasLineElement(SUSPECT, true));
        assertTrue(bnd2.hasLineElement(SUSPECT, true));
        assertFalse(bnd1.isLineElement(false));
        assertFalse(bnd2.isLineElement(false));

        cluster = tester.Analyser.getClusters().get(0);
        assertTrue(cluster.hasLinkingLineElements());

        // now test BNDs in a DB on the LINE arm which will invalidate the line test
        bnd1 = createSv(tester.nextVarId(), "1", "2", 100, 100,  1, 1, BND, "");
        bnd2 = createSv(tester.nextVarId(), "1", "2", 110, 110,  -1, -1, BND, polyTMotif);

        tester.addAndCluster(bnd1, bnd2);

        assertFalse(bnd1.isLineElement(true));
        assertFalse(bnd2.isLineElement(true));
        assertFalse(bnd1.isLineElement(false));
        assertFalse(bnd2.isLineElement(false));

        cluster = tester.Analyser.getClusters().get(0);
        assertTrue(!cluster.hasLinkingLineElements());

        // some other variant having the poly A motif
        tester.clearClustersAndSVs();

        bnd1 = createSv(tester.nextVarId(), "1", "2", 100, 100,  -1, 1, BND, "");
        bnd2 = createSv(tester.nextVarId(), "1", "3", 150, 200,  -1, -1, BND, "");
        SvVarData del = createSv(tester.nextVarId(), "1", "1", 1000, 11000,  1, -1, DEL, polyAMotif);

        tester.AllVariants.addAll(Lists.newArrayList(bnd1, bnd2, del));
        tester.preClusteringInit();
        tester.Analyser.clusterAndAnalyse();

        assertTrue(bnd1.hasLineElement(SUSPECT, true));
        assertTrue(bnd2.hasLineElement(SUSPECT, true));
        assertTrue(del.hasLineElement(SUSPECT, true));

        assertFalse(del.isLineElement(false));
        assertFalse(bnd1.isLineElement(false));
        assertFalse(bnd2.isLineElement(false));

        cluster = tester.Analyser.getClusters().get(0);
        assertTrue(cluster.hasLinkingLineElements());

        tester.clearClustersAndSVs();

        // more than 1 BND and a SGL
        String insertPolyTMotif = shortBases + POLY_T_MOTIF + randomBases;
        sgl = createSv(tester.nextVarId(), "1", "", 100, 0,  1, 0, SGL, insertPolyTMotif);
        bnd1 = createSv(tester.nextVarId(), "1", "2", 110, 100,  -1, 1, BND, "");
        bnd2 = createSv(tester.nextVarId(), "2", "3", 200, 1000,  1, -1, BND, "");

        tester.AllVariants.addAll(Lists.newArrayList(bnd1, bnd2, sgl));
        tester.preClusteringInit();
        tester.Analyser.clusterAndAnalyse();

        assertTrue(bnd1.hasLineElement(SUSPECT, false));
        assertTrue(bnd2.hasLineElement(SUSPECT, true));
        assertTrue(!sgl.inLineElement());

        assertFalse(bnd1.isLineElement(true));
        assertFalse(bnd2.isLineElement(false));

        cluster = tester.Analyser.getClusters().get(0);
        assertTrue(cluster.hasLinkingLineElements());

        // BNDs which would be suspect except for the remote DB falling in a known line location
        bnd1 = createSv(tester.nextVarId(), "1", "2", 100, 100,  -1, 1, BND, "");
        bnd2 = createSv(tester.nextVarId(), "1", "2", 200, 110,  1, -1, BND, polyAMotif);
        bnd1.addLineElement(KNOWN, false);

        tester.addAndCluster(bnd1, bnd2);

        assertFalse(bnd1.hasLineElement(SUSPECT, true));
        assertFalse(bnd2.hasLineElement(SUSPECT, true));
        assertFalse(bnd1.isLineElement(true));
        assertFalse(bnd2.isLineElement(false));

        cluster = tester.Analyser.getClusters().get(0);
        assertTrue(!cluster.hasLinkingLineElements());

        // long INV plus a single in a remote DB
        SvVarData inv = createSv(tester.nextVarId(), "1", "1", 100, 1050000,  1, 1, INV, polyAMotif);
        sgl = createSv(tester.nextVarId(), "1", "0", 1050020, 0,  -1, 0, SGL, polyAMotif);

        tester.addAndCluster(inv, sgl);

        assertTrue(inv.hasLineElement(SUSPECT, true));
        assertFalse(inv.hasLineElement(SUSPECT, false));

        cluster = tester.Analyser.getClusters().get(0);
        assertTrue(cluster.hasLinkingLineElements());
    }

    @Test
    public void testSuspectSglLineMarking()
    {
        LinxTester tester = new LinxTester();

        String shortBases = "GCTG";
        String randomBases = "GCTGTGAGCTGATCG";
        String insertPolyAMotif = randomBases + POLY_A_MOTIF + shortBases;

        // 2 SGLs with an insert-site motif and in a short DB
        SvVarData sgl1 = createSv(tester.nextVarId(), "1", "", 100, 0,  1, 0, SGL, "");
        SvVarData sgl2 = createSv(tester.nextVarId(), "1", "", 120, 0,  -1, 0, SGL, insertPolyAMotif);

        tester.addAndCluster(sgl1, sgl2);

        assertFalse(sgl1.isLineElement(false));
        assertFalse(sgl2.isLineElement(false));

        SvCluster cluster = tester.Analyser.getClusters().get(0);
        assertTrue(cluster.hasLinkingLineElements());

        tester.clearClustersAndSVs();

        // can also be just with a lone SGL as long as the CN change is relatively small
        sgl1 = createSv(tester.nextVarId(), "1", "", 120, 0,  -1, 0, SGL, insertPolyAMotif);

        tester.AllVariants.add(sgl1);
        tester.preClusteringInit();
        tester.Analyser.clusterAndAnalyse();

        assertFalse(sgl1.isLineElement(false));

        cluster = tester.Analyser.getClusters().get(0);
        assertTrue(!cluster.hasLinkingLineElements());

        tester.clearClustersAndSVs();

        // now with a low CN change
        sgl1 = createTestSv(tester.nextVarId(), "1", "", 120, 0,  -1, 0, SGL,
                2, 0, 0.1, 0, 0.1, insertPolyAMotif);

        tester.AllVariants.add(sgl1);
        tester.preClusteringInit();
        tester.Analyser.clusterAndAnalyse();

        assertFalse(sgl1.isLineElement(false));

        cluster = tester.Analyser.getClusters().get(0);
        assertTrue(cluster.hasLinkingLineElements());
    }

    @Test
    public void testLineChaining()
    {
        LinxTester tester = new LinxTester();

        // Configurator.setRootLevel(Level.DEBUG);

        // scenario 1: simple line chain with assembled breakends
        SvVarData var1 = createBnd(tester.nextVarId(), "1", 100, -1, "2", 100, 1);
        SvVarData var2 = createBnd(tester.nextVarId(), "1", 200, 1, "2", 200, -1);
        SvVarData var3 = createBnd(tester.nextVarId(), "1", 1100, -1, "3", 100, 1);
        SvVarData var4 = createBnd(tester.nextVarId(), "1", 1200, 1, "3", 200, -1);

        var1.addLineElement(KNOWN, true);
        var2.addLineElement(KNOWN, true);
        var3.addLineElement(KNOWN, true);
        var4.addLineElement(KNOWN, true);
        var1.setAssemblyData(true, "asmb12");
        var2.setAssemblyData(true, "asmb12");

        tester.AllVariants.addAll(Lists.newArrayList(var1, var2, var3, var4));

        tester.preClusteringInit();
        tester.Analyser.clusterAndAnalyse();

        assertEquals(1, tester.Analyser.getClusters().size());
        SvCluster cluster = tester.Analyser.getClusters().get(0);

        assertEquals(2, cluster.getChains().size());

        tester.clearClustersAndSVs();

        // more complicated scenario involving multiple source locations and breakends at each location

        var1 = createBnd(tester.nextVarId(), "1", 100, -1, "3", 100, 1);
        var2 = createBnd(tester.nextVarId(), "1", 200, 1, "3", 200, -1);

        var3 = createBnd(tester.nextVarId(), "1", 300, -1, "2", 100, 1);
        var4 = createInv(tester.nextVarId(), "1", 400, 600, 1);
        SvVarData var5 = createBnd(tester.nextVarId(), "1", 500, -1, "2", 200, -1);

        var3.setAssemblyData(true, "asmb12");
        var4.setAssemblyData(true, "asmb12");

        // another independent source element
        SvVarData var6 = createBnd(tester.nextVarId(), "3", 1000, -1, "4", 100, 1);
        SvVarData var7 = createBnd(tester.nextVarId(), "3", 1100, 1, "4", 200, -1);

        // and another but with breakends going to different insert locations, which is invalid
        SvVarData var8 = createBnd(tester.nextVarId(), "1", 2000, -1, "5", 100, 1);
        SvVarData var9 = createBnd(tester.nextVarId(), "1", 2100, 1, "6", 200, -1);

        tester.AllVariants.addAll(Lists.newArrayList(var1, var2, var3, var4,  var5, var6, var7, var8, var9));

        tester.AllVariants.forEach(x -> x.addLineElement(KNOWN, true));

        tester.preClusteringInit();
        tester.Analyser.clusterAndAnalyse();

        assertEquals(1, tester.Analyser.getClusters().size());
        cluster = tester.Analyser.getClusters().get(0);

        assertEquals(3, cluster.getChains().size());
        assertTrue(tester.findChainWithSVs(cluster, Lists.newArrayList(var1, var2)) != null);
        assertTrue(tester.findChainWithSVs(cluster, Lists.newArrayList(var3, var4, var5)) != null);
        assertTrue(tester.findChainWithSVs(cluster, Lists.newArrayList(var6, var7)) != null);
    }

    @Test
    public void testLineSglChaining()
    {
        LinxTester tester = new LinxTester();

        // Configurator.setRootLevel(Level.DEBUG);

        String shortBases = "GCTG";
        String randomBases = "GCTGTGAGCTGATCG";
        String insertPolyTMotif = shortBases + POLY_T_MOTIF + randomBases;

        // scenario 1: 2 SGLs with mappings to the same location
        SvVarData sgl1 = createSv(tester.nextVarId(), "1", "", 100, 0, POS_ORIENT, 0, SGL, insertPolyTMotif);
        SvVarData sgl2 = createSv(tester.nextVarId(), "1", "", 120, 0, NEG_ORIENT, 0, SGL, "");

        // both have additional mappings which are ignored
        sgl1.getSglMappings().add(new SglMapping(CHR_2, 1000, NEG_ORIENT, "", 1));
        sgl1.getSglMappings().add(new SglMapping("5", 10000, NEG_ORIENT, "", 1));
        sgl2.getSglMappings().add(new SglMapping(CHR_2, 1100, POS_ORIENT, "", 1));
        sgl2.getSglMappings().add(new SglMapping("3", 2000, NEG_ORIENT, "", 1));

        tester.AllVariants.addAll(Lists.newArrayList(sgl1, sgl2));

        tester.preClusteringInit();
        tester.Analyser.clusterAndAnalyse();

        assertEquals(1, tester.Analyser.getClusters().size());
        SvCluster cluster = tester.Analyser.getClusters().get(0);

        assertEquals(1, cluster.getChains().size());

        tester.clearClustersAndSVs();

        // again with a BND to a SGL
        sgl1 = createSv(tester.nextVarId(), "1", "", 100, 0,  1, 0, SGL, insertPolyTMotif);
        SvVarData bnd = createSv(tester.nextVarId(), "1", "2", 120, 1100,  -1, 1, BND, "");
        bnd.addLineElement(KNOWN, false);

        // both have additional mappings which are ignored
        sgl1.getSglMappings().add(new SglMapping(CHR_2, 1000, NEG_ORIENT, "", 1));
        sgl1.getSglMappings().add(new SglMapping("5", 10000, NEG_ORIENT, "", 1));

        tester.AllVariants.addAll(Lists.newArrayList(sgl1, bnd));

        tester.preClusteringInit();
        tester.Analyser.clusterAndAnalyse();

        assertEquals(1, tester.Analyser.getClusters().size());
        cluster = tester.Analyser.getClusters().get(0);
        assertTrue(cluster.hasLinkingLineElements());

        assertEquals(1, cluster.getChains().size());
    }

}