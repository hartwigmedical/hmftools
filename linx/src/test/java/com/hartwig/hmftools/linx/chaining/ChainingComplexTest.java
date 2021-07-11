package com.hartwig.hmftools.linx.chaining;

import static com.hartwig.hmftools.common.sv.StructuralVariantType.BND;
import static com.hartwig.hmftools.common.sv.StructuralVariantType.DEL;
import static com.hartwig.hmftools.common.sv.StructuralVariantType.DUP;
import static com.hartwig.hmftools.common.sv.StructuralVariantType.INV;
import static com.hartwig.hmftools.linx.utils.SvTestUtils.createTestSv;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertTrue;

import com.hartwig.hmftools.linx.types.SvCluster;
import com.hartwig.hmftools.linx.types.SvVarData;
import com.hartwig.hmftools.linx.utils.LinxTester;

import org.junit.Ignore;
import org.junit.Test;

public class ChainingComplexTest
{
    @Test
    public void testBFBChain1()
    {
        // vanilla BFB of the form centromere - 1 - 2 - 1 - 3 - 1 - 2 - 1 - 4 - 1 - 2 - 1 - 3 - 1 - 2 - 1 - R - telomere, where R is the resolving SV
        LinxTester tester = new LinxTester();

        // Configurator.setRootLevel(Level.DEBUG);

        final SvVarData var1 = createTestSv(1, "1", "1", 1000,2000, -1, -1, INV,  8);
        final SvVarData var2 = createTestSv(2, "1", "1", 9000,10000, 1, 1, INV,  4);
        final SvVarData var3 = createTestSv(3, "1", "1", 6000,7000, 1, 1, INV, 2);
        final SvVarData var4 = createTestSv(4, "1", "1", 3000,4000, 1, 1, INV, 1);
        final SvVarData var5 = createTestSv(5, "1", "2", 12000,100, 1, 1, BND, 1);

        tester.AllVariants.add(var1);
        tester.AllVariants.add(var2);
        tester.AllVariants.add(var3);
        tester.AllVariants.add(var4);
        tester.AllVariants.add(var5);

        tester.preClusteringInit();
        tester.Analyser.clusterAndAnalyse();

        assertEquals(var1.getFoldbackId(true), var1.id());
        assertEquals(var2.getFoldbackId(true), var2.id());
        assertEquals(var3.getFoldbackId(true), var3.id());
        assertEquals(var4.getFoldbackId(true), var4.id());

        assertEquals(1, tester.Analyser.getClusters().size());
        final SvCluster cluster = tester.Analyser.getClusters().get(0);

        assertEquals(1, cluster.getChains().size());

        final SvChain chain = cluster.getChains().get(0);

        assertEquals(5, chain.getSvCount());
        assertEquals(15, chain.getLinkCount());
    }

    @Test
    public void testBFBChain2()
    {
        // vanilla BFB of the form centromere - 1 - 2 - 1 - 3 - 1 - 4 - telomere, where R is the resolving SV
        LinxTester tester = new LinxTester();

        final SvVarData var1 = createTestSv(1, "1", "1", 2000,3000, -1, -1, INV,  3);
        final SvVarData var2 = createTestSv(2, "1", "1", 9000,10000, 1, 1, INV,  1);
        final SvVarData var3 = createTestSv(3, "1", "1", 5000,6000, 1, 1, INV, 1);
        final SvVarData var4 = createTestSv(4, "1", "1", 1000,8000, 1, 1, INV, 1);

        tester.AllVariants.add(var1);
        tester.AllVariants.add(var2);
        tester.AllVariants.add(var3);
        tester.AllVariants.add(var4);

        tester.preClusteringInit();
        tester.Analyser.clusterAndAnalyse();

        assertEquals(var1.getFoldbackId(true), var1.id());
        assertEquals(var2.getFoldbackId(true), var2.id());
        assertEquals(var3.getFoldbackId(true), var3.id());

        assertEquals(1, tester.Analyser.getClusters().size());
        final SvCluster cluster = tester.Analyser.getClusters().get(0);

        assertEquals(1, cluster.getChains().size());

        final SvChain chain = cluster.getChains().get(0);

        assertEquals(5, chain.getLinkCount());
    }

    @Test
    public void testBFBChainWithChainedFoldbacks()
    {
        // BFB of the form centromere - A01 - B2 - A01 - C3 - A01 - B2 - A01 - R4 - telomere, where R is the resolving SV
        // but assembled fragments are inserted into each of the foldbacks, requiring linking to make use of matching ploidy
        LinxTester tester = new LinxTester();

        final SvVarData varA1 = createTestSv(0, "1", "1", 2000,5500, -1, -1, INV, 4);
        final SvVarData varA2 = createTestSv(1, "1", "1", 3000,5600, -1, 1, DUP, 4);

        varA1.setAssemblyData(false, "asmb_A1_A2");
        varA2.setAssemblyData(false, "asmb_A1_A2");

        final SvVarData varB = createTestSv(2, "1", "1", 9000,10000, 1, 1, INV, 2);

        // functions as a foldback but interrupted by the assembled TI A1-A2
        final SvVarData varC = createTestSv(3, "1", "1", 5000,6000, 1, 1, INV, 1);
        final SvVarData varR = createTestSv(4, "1", "1", 1000,8000, 1, 1, INV, 1);

        // CN profile
        // T - 2 - 1 - 5 - 9 - 8 - 12 - 8 - 7 - 6 - 4 - 2 - C

        tester.AllVariants.add(varA1);
        tester.AllVariants.add(varA2);
        tester.AllVariants.add(varB);
        tester.AllVariants.add(varC);
        tester.AllVariants.add(varR);

        tester.preClusteringInit();
        tester.Analyser.clusterAndAnalyse();

        assertTrue(varB.isFoldback());
        assertTrue(varC.isFoldback());
        assertTrue(varA1.isChainedFoldback());
        assertTrue(varA2.isChainedFoldback());
        assertEquals(varA2.id(), varA1.getFoldbackId(true));
        assertEquals(varA1.id(), varA2.getFoldbackId(true));
        assertEquals(varB.id(), varB.getFoldbackId(true));

        assertEquals(1, tester.Analyser.getClusters().size());
        final SvCluster cluster = tester.Analyser.getClusters().get(0);

        assertEquals(1, cluster.getChains().size());

        final SvChain chain = cluster.getChains().get(0);

        int linkCount = 2 * 4 + 1 * 2 + 2 - 1;
        assertEquals(linkCount, chain.getLinkCount());
    }

    @Test
    public void testComplexDupAssemblyChain()
    {
        // simple chain with replicated section: telomere - DEL - remote TI - DUP - DEL - remote TI - centromere
        // SVs: telomere - 1 -
        LinxTester tester = new LinxTester();

        final SvVarData var1 = createTestSv(1, "1", "1", 2000,3000, 1, -1, DEL, 2);
        final SvVarData var2 = createTestSv(2, "1", "2", 4000,100, 1, -1, BND, 2);
        final SvVarData var3 = createTestSv(3, "1", "2", 5000,200, -1, 1, BND, 2);
        final SvVarData dup = createTestSv(4, "1", "1", 1000,6000, -1, 1, DUP, 1);

        var1.setAssemblyData(false, "asmb12");
        var2.setAssemblyData(true, "asmb12");

        var2.setAssemblyData(false, "asmb23");
        var3.setAssemblyData(false, "asmb23");

        tester.AllVariants.add(var1);
        tester.AllVariants.add(var2);
        tester.AllVariants.add(var3);
        tester.AllVariants.add(dup);

        tester.preClusteringInit(true);
        tester.Analyser.clusterAndAnalyse();

        assertEquals(1, tester.Analyser.getClusters().size());
        final SvCluster cluster = tester.Analyser.getClusters().get(0);

        assertTrue(cluster.getFoldbacks().isEmpty());
        assertEquals(1, cluster.getChains().size());

        final SvChain chain = cluster.getChains().get(0);

        assertEquals(6, chain.getLinkCount());
    }

    @Test
    public void testComplexMultiAssemblyChain()
    {
        // chain with more than breakend connecting via assembly to other breakends
        LinxTester tester = new LinxTester();

        final SvVarData var1 = createTestSv(1, "1", "1", 2000,1000, 1, -1, DEL, 5);
        final SvVarData var2 = createTestSv(2, "1", "1", 1500,9000, 1, 1, INV, 5);
        final SvVarData var3 = createTestSv(3, "1", "1", 1200,6000, 1, -1, DEL, 4);
        final SvVarData var4 = createTestSv(4, "1", "2", 1300,200, -1, 1, BND, 2);
        final SvVarData var5 = createTestSv(5, "1", "2", 1400,500, -1, -1, BND, 1);
        final SvVarData var6 = createTestSv(6, "1", "1", 100,900, 1, -1, DEL, 2);
        final SvVarData var7 = createTestSv(7, "1", "1", 800,2000, -1, -1, INV, 1);
        final SvVarData var8 = createTestSv(8, "1", "1", 1100,10000, 1, 1, INV, 2);

        var1.setAssemblyData(false, "asmb12;asmb18;asmb13");
        var2.setAssemblyData(true, "asmb12;asmb24;asmb25");
        var3.setAssemblyData(true, "asmb13;asmb36;asmb37");
        var4.setAssemblyData(true, "asmb24");
        var5.setAssemblyData(true, "asmb25");
        var6.setAssemblyData(false, "asmb36");
        var7.setAssemblyData(true, "asmb37");
        var8.setAssemblyData(true, "asmb18");

        tester.AllVariants.add(var1);
        tester.AllVariants.add(var2);
        tester.AllVariants.add(var3);
        tester.AllVariants.add(var4);
        tester.AllVariants.add(var5);
        tester.AllVariants.add(var6);
        tester.AllVariants.add(var7);
        tester.AllVariants.add(var8);

        tester.preClusteringInit(true);
        tester.Analyser.clusterAndAnalyse();

        assertEquals(1, tester.Analyser.getClusters().size());
        final SvCluster cluster = tester.Analyser.getClusters().get(0);

        assertEquals(7, cluster.getAssemblyLinkedPairs().size());

        int assemblyChainedLinks = cluster.getChains().stream().mapToInt(x -> x.getAssemblyLinkCount()).sum();
        assertEquals(7, assemblyChainedLinks);
    }

    @Test
    public void testBFBChainWith2ChainedFoldbacks()
    {
        // BFB of the form telomere - 3 - 4 - 3 - 1 - 2 - 1 - 3 - 4 - 3 - R5-6 - centromere, where R is the resolving SV
        LinxTester tester = new LinxTester();

        final SvVarData var1 = createTestSv(1, "1", "1", 1000,2000, -1, -1, INV, 2);
        final SvVarData var2 = createTestSv(2, "1", "1", 3000,4000, 1, 1, INV, 1);
        final SvVarData var3 = createTestSv(3, "1", "1", 6000,7000, 1, -1, DEL, 4);
        final SvVarData var4 = createTestSv(4, "1", "1", 9000,10000, 1, 1, INV, 2);
        final SvVarData var5 = createTestSv(5, "1", "1", 5000,15000, -1, -1, INV, 1);
        final SvVarData var6 = createTestSv(6, "1", "2", 15500,1000, 1, 1, BND, 1);

        var1.setAssemblyData(true, "asmb12s;asmb12e");
        var2.setAssemblyData(true, "asmb12s");
        var2.setAssemblyData(false, "asmb12e");

        var3.setAssemblyData(false, "asmb34s;asmb34e");
        var4.setAssemblyData(true, "asmb34s");
        var4.setAssemblyData(false, "asmb34e");

        tester.AllVariants.add(var1);
        tester.AllVariants.add(var2);
        tester.AllVariants.add(var3);
        tester.AllVariants.add(var4);
        tester.AllVariants.add(var5);
        tester.AllVariants.add(var6);

        tester.preClusteringInit();
        tester.Analyser.clusterAndAnalyse();

        assertTrue(var1.isFoldback());
        assertTrue(var3.isFoldback());

        assertEquals(1, tester.Analyser.getClusters().size());
        final SvCluster cluster = tester.Analyser.getClusters().get(0);

        assertEquals(1, cluster.getChains().size());

        final SvChain chain = cluster.getChains().get(0);

        assertEquals(10, chain.getLinkCount());
        assertEquals(6, chain.getSvCount());
        assertEquals(1, chain.jcn(), 0.01);
    }

    @Test
    public void testBFBChainToComplexChain()
    {
        // chained foldback connects to higher ploidy chained section
        //  T - 3 - 4 - 1 - 2 - 3 - 4 - 5 - C
        LinxTester tester = new LinxTester();

        final SvVarData var1 = createTestSv(1, "1", "2", 1000,1000, -1, -1, BND, 1);
        final SvVarData var2 = createTestSv(2, "1", "2", 6000,3000, -1, 1, BND, 1);

        final SvVarData var3 = createTestSv(3, "1", "2", 10000,20000, 1, -1, BND, 2);
        final SvVarData var4 = createTestSv(4, "1", "2", 11000,25000, -1, 1, BND, 2);
        final SvVarData var5 = createTestSv(5, "1", "3", 28000,2000, 1, 1, BND, 1);

        // mark breakends as assembled to ensure they are foldbacks
        var1.setAssemblyData(false, "asmb12");
        var2.setAssemblyData(false, "asmb12");

        var3.setAssemblyData(false, "asmb34");
        var4.setAssemblyData(false, "asmb34");

        tester.AllVariants.add(var1);
        tester.AllVariants.add(var2);
        tester.AllVariants.add(var3);
        tester.AllVariants.add(var4);
        tester.AllVariants.add(var5);

        tester.preClusteringInit(true);
        tester.Analyser.clusterAndAnalyse();

        assertTrue(var1.isFoldback());
        assertTrue(var2.isFoldback());

        assertEquals(1, tester.Analyser.getClusters().size());
        final SvCluster cluster = tester.Analyser.getClusters().get(0);

        assertEquals(1, cluster.getChains().size());

        final SvChain chain = cluster.getChains().get(0);

        assertEquals(6, chain.getLinkCount());
        assertEquals(5, chain.getSvCount());
        assertEquals(1, chain.jcn(), 0.01);
    }

    @Test
    public void testComplexDupSimple()
    {
        // simple chain with replicated SV: telomere - 4 - 5 - 1 - 2 - DUP - 1 - 3 - centromere
        LinxTester tester = new LinxTester();

        final SvVarData var1 = createTestSv(1, "1", "1", 4000,5000, 1, -1, DEL, 2);
        final SvVarData var2 = createTestSv(2, "1", "1", 2000,6000, -1, 1, DUP, 1);
        final SvVarData var3 = createTestSv(3, "1", "1", 7000,8000, 1, -1, DEL, 1);

        // add some BND so the group isn't considered simple and split up
        final SvVarData var4 = createTestSv(4, "1", "2", 500,100, 1, -1, BND, 1);
        final SvVarData var5 = createTestSv(5, "1", "2", 1000,200, -1, 1, BND, 1);

        tester.AllVariants.add(var1);
        tester.AllVariants.add(var2);
        tester.AllVariants.add(var3);
        tester.AllVariants.add(var4);
        tester.AllVariants.add(var5);

        tester.preClusteringInit();
        tester.Analyser.clusterAndAnalyse();

        assertEquals(1, tester.Analyser.getClusters().size());
        final SvCluster cluster = tester.Analyser.getClusters().get(0);

        assertTrue(cluster.getFoldbacks().isEmpty());
        assertEquals(1, cluster.getChains().size());

        final SvChain chain = cluster.getChains().get(0);

        assertEquals(5, chain.getLinkCount());
        assertEquals(5, chain.getSvCount());
    }

    @Test
    public void testComplexDupChain()
    {
        // simple chain with replicated section: centromere - A - B - C - D - DUP - B - C - D - E - telomere,
        // where D is a complex DUP around the section B - C - D
        LinxTester tester = new LinxTester();

        final SvVarData varA = createTestSv(1, "1", "1", 1000,5000, 1, 1, INV, 1);
        final SvVarData varB = createTestSv(2, "1", "1", 4000,9000, -1, 1, DUP, 2);
        final SvVarData varC = createTestSv(3, "1", "1", 7000,8000, 1, -1, DEL, 2);
        final SvVarData varD = createTestSv(4, "1", "1", 3000,6000, 1, -1, DEL, 2);
        final SvVarData varE = createTestSv(5, "1", "1", 2000,10000, -1, -1, INV, 1);
        final SvVarData varDup = createTestSv(6, "1", "1", 2500,4500, -1, 1, DUP, 1);

        tester.AllVariants.add(varA);
        tester.AllVariants.add(varB);
        tester.AllVariants.add(varC);
        tester.AllVariants.add(varD);
        tester.AllVariants.add(varE);
        tester.AllVariants.add(varDup);

        tester.preClusteringInit();
        tester.Analyser.clusterAndAnalyse();

        assertEquals(1, tester.Analyser.getClusters().size());
        final SvCluster cluster = tester.Analyser.getClusters().get(0);

        assertTrue(cluster.getFoldbacks().isEmpty());
        assertEquals(1, cluster.getChains().size());

        final SvChain chain = cluster.getChains().get(0);

        assertEquals(8, chain.getLinkCount());
    }

    @Ignore
    @Test
    // difficult for chaining to make the correct links
    public void testBFBChainWithComplexDup()
    {
        // BFB of the form centromere - 1 - 2 - 1 - 3 - 1 - 4-DUP - 2 - 1 - 3 - 1 - 2 - 5R - telomere,
        // where D is a complex DUP around the section B - A - C - A and R is the resolving SV
        LinxTester tester = new LinxTester();

        final SvVarData var1 = createTestSv(1, "1", "1", 2000,3000, -1, -1, INV, 5);
        final SvVarData var2 = createTestSv(2, "1", "1", 9000,10000, 1, 1, INV, 3);
        final SvVarData var3 = createTestSv(3, "1", "1", 5000,6000, 1, 1, INV, 2);
        final SvVarData var4 = createTestSv(4, "1", "1", 7000,8000, -1, 1, DUP, 1);
        final SvVarData var5 = createTestSv(5, "1", "1", 1000,4000, 1, -1, DEL, 1);

        tester.AllVariants.add(var1);
        tester.AllVariants.add(var2);
        tester.AllVariants.add(var3);
        tester.AllVariants.add(var4);
        tester.AllVariants.add(var5);

        tester.preClusteringInit();
        tester.Analyser.clusterAndAnalyse();

        assertEquals(var1.getFoldbackId(true), var1.id());
        assertEquals(var2.getFoldbackId(true), var2.id());
        assertEquals(var3.getFoldbackId(true), var3.id());

        assertEquals(1, tester.Analyser.getClusters().size());
        final SvCluster cluster = tester.Analyser.getClusters().get(0);

        assertEquals(1, cluster.getChains().size());

        final SvChain chain = cluster.getChains().get(0);

        assertEquals(11, chain.getLinkCount());
    }

}
