package com.hartwig.hmftools.linx.analyser;

import static com.hartwig.hmftools.common.variant.structural.StructuralVariantType.BND;
import static com.hartwig.hmftools.common.variant.structural.StructuralVariantType.DEL;
import static com.hartwig.hmftools.common.variant.structural.StructuralVariantType.DUP;
import static com.hartwig.hmftools.common.variant.structural.StructuralVariantType.INV;
import static com.hartwig.hmftools.linx.analyser.SvTestHelper.createDel;
import static com.hartwig.hmftools.linx.analyser.SvTestHelper.createDup;
import static com.hartwig.hmftools.linx.analyser.SvTestHelper.createInv;
import static com.hartwig.hmftools.linx.analyser.SvTestHelper.createTestSv;
import static com.hartwig.hmftools.linx.chaining.ChainFinder.CHAIN_METHOD_COMPARE;
import static com.hartwig.hmftools.linx.chaining.ChainFinder.CHAIN_METHOD_NEW;
import static com.hartwig.hmftools.linx.chaining.ChainFinder.CHAIN_METHOD_OLD;
import static com.hartwig.hmftools.linx.chaining.SvChain.CHAIN_ASSEMBLY_LINK_COUNT;
import static com.hartwig.hmftools.linx.chaining.SvChain.CHAIN_LINK_COUNT;
import static com.hartwig.hmftools.linx.chaining.SvChain.checkIsValid;
import static com.hartwig.hmftools.linx.types.SvLinkedPair.LINK_TYPE_TI;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertFalse;
import static org.junit.Assert.assertTrue;

import com.hartwig.hmftools.linx.chaining.SvChain;
import com.hartwig.hmftools.linx.types.SvCluster;
import com.hartwig.hmftools.linx.types.SvLinkedPair;
import com.hartwig.hmftools.linx.types.SvVarData;

import org.junit.Ignore;
import org.junit.Test;

public class ChainingComplexTest
{

    @Test
    public void testBFBChain1()
    {
        // vanilla BFB of the form centromere - A - B - A - C - A - R - telomere, where R is the resolving SV
        SvTestHelper tester = new SvTestHelper();
        tester.logVerbose(true);

        final SvVarData varA = createTestSv("0", "1", "1", 2000,3000, -1, -1, INV,  3);
        final SvVarData varB = createTestSv("1", "1", "1", 9000,10000, 1, 1, INV,  1);
        final SvVarData varC = createTestSv("2", "1", "1", 5000,6000, 1, 1, INV, 1);
        final SvVarData varR = createTestSv("3", "1", "1", 1000,8000, 1, 1, INV, 1);

        tester.AllVariants.add(varA);
        tester.AllVariants.add(varB);
        tester.AllVariants.add(varC);
        tester.AllVariants.add(varR);

        tester.preClusteringInit();
        tester.Analyser.clusterAndAnalyse();

        assertEquals(varA.getFoldbackLink(true), varA.id());
        assertEquals(varB.getFoldbackLink(true), varB.id());
        assertEquals(varC.getFoldbackLink(true), varC.id());

        assertEquals(1, tester.Analyser.getClusters().size());
        final SvCluster cluster = tester.Analyser.getClusters().get(0);

        assertEquals(1, cluster.getChains().size());

        final SvChain chain = cluster.getChains().get(0);

        assertEquals(5, chain.getLinkCount());
    }

    @Ignore
    @Test
    public void testBFBChainWithComplexDup()
    {
        // BFB of the form centromere - 1 - 2 - 1 - 3 - 1 - 4-DUP - 2 - 1 - 3 - 1 - 2 - 5R - telomere,
        // where D is a complex DUP around the section B - A - C - A and R is the resolving SV
        SvTestHelper tester = new SvTestHelper();
        tester.logVerbose(true);

        // tester.setChaininMethod(CHAIN_METHOD_OLD);

        final SvVarData var1 = createTestSv("1", "1", "1", 2000,3000, -1, -1, INV, 5);
        final SvVarData var2 = createTestSv("2", "1", "1", 9000,10000, 1, 1, INV, 3);
        final SvVarData var3 = createTestSv("3", "1", "1", 5000,6000, 1, 1, INV, 2);
        final SvVarData var4 = createTestSv("4", "1", "1", 7000,8000, -1, 1, DUP, 1);
        final SvVarData var5 = createTestSv("5", "1", "1", 1000,4000, 1, -1, DEL, 1);

        // CN profile:
        // T - 2 - 1 - 6 - 11 - 12 - 10 - 8 - 9 - 8 - 5 - 2 - C

        tester.AllVariants.add(var1);
        tester.AllVariants.add(var2);
        tester.AllVariants.add(var3);
        tester.AllVariants.add(var4);
        tester.AllVariants.add(var5);

        tester.preClusteringInit();
        tester.Analyser.clusterAndAnalyse();

        assertEquals(var1.getFoldbackLink(true), var1.id());
        assertEquals(var2.getFoldbackLink(true), var2.id());
        assertEquals(var3.getFoldbackLink(true), var3.id());

        assertEquals(1, tester.Analyser.getClusters().size());
        final SvCluster cluster = tester.Analyser.getClusters().get(0);

        assertEquals(1, cluster.getChains().size());

        final SvChain chain = cluster.getChains().get(0);

        assertEquals(11, chain.getLinkCount());
    }

    @Test
    public void testComplexDupSimple()
    {
        // simple chain with replicated SV: telomere - 4 - 5 - 1 - 2-DUP - 1 - 3 - centromere
        SvTestHelper tester = new SvTestHelper();
        tester.logVerbose(true);

        final SvVarData var1 = createTestSv("1", "1", "1", 4000,5000, 1, -1, DEL, 2);
        final SvVarData var2 = createTestSv("2", "1", "1", 2000,6000, -1, 1, DUP, 1);
        final SvVarData var3 = createTestSv("3", "1", "1", 7000,8000, 1, -1, DEL, 1);

        // add some BND so the group isn't considered simple and split up
        final SvVarData var4 = createTestSv("4", "1", "2", 500,100, 1, -1, BND, 1);
        final SvVarData var5 = createTestSv("5", "1", "2", 1000,200, -1, 1, BND, 1);

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
    }

    @Test
    public void testComplexDupChain()
    {
        // simple chain with replicated section: centromere - A - B - C - D - DUP - B - C - D - E - telomere,
        // where D is a complex DUP around the section B - C - D
        SvTestHelper tester = new SvTestHelper();
        tester.logVerbose(false);

        final SvVarData varA = createTestSv("1", "1", "1", 1000,5000, 1, 1, INV, 1);
        final SvVarData varB = createTestSv("2", "1", "1", 4000,9000, -1, 1, DUP, 2);
        final SvVarData varC = createTestSv("3", "1", "1", 7000,8000, 1, -1, DEL, 2);
        final SvVarData varD = createTestSv("4", "1", "1", 3000,6000, 1, -1, DEL, 2);
        final SvVarData varE = createTestSv("5", "1", "1", 2000,10000, -1, -1, INV, 1);
        final SvVarData varDup = createTestSv("6", "1", "1", 2500,4500, -1, 1, DUP, 1);

        // CN profile:
        // T - 2 - 1 - 2 - 3 - 1 - 3 - 2 - 1 - 3 - 1 - 3 - 1 - 2  - C

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

    @Test
    public void testBFBChainWithChainedFoldbacks()
    {
        // BFB of the form centromere - A01 - B2 - A01 - C3 - A01 - B2 - A01 - R4 - telomere, where R is the resolving SV
        // but assembled fragments are inserted into each of the foldbacks, requiring linking to make use of matching ploidy
        SvTestHelper tester = new SvTestHelper();

        // tester.setChaininMethod(CHAIN_METHOD_COMPARE);
        tester.logVerbose(true);

        final SvVarData varA1 = createTestSv("0", "1", "1", 2000,5500, -1, -1, INV, 4);
        final SvVarData varA2 = createTestSv("1", "1", "1", 3000,5600, -1, 1, DUP, 4);

        varA1.setAssemblyData(false, "asmb_A1_A2");
        varA2.setAssemblyData(false, "asmb_A1_A2");

        final SvVarData varB = createTestSv("2", "1", "1", 9000,10000, 1, 1, INV, 2);

        // functions as a foldback but interrupted by the assembled TI A1-A2
        final SvVarData varC = createTestSv("3", "1", "1", 5000,6000, 1, 1, INV, 1);
        final SvVarData varR = createTestSv("4", "1", "1", 1000,8000, 1, 1, INV, 1);

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
        assertEquals(varA2.id(), varA1.getFoldbackLink(true));
        assertEquals(varA1.id(), varA2.getFoldbackLink(true));
        assertEquals(varB.id(), varB.getFoldbackLink(true));

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
        SvTestHelper tester = new SvTestHelper();
        tester.logVerbose(true);

        final SvVarData var1 = createTestSv("1", "1", "1", 2000,3000, 1, -1, DEL, 2);
        final SvVarData var2 = createTestSv("2", "1", "2", 4000,100, 1, -1, BND, 2);
        final SvVarData var3 = createTestSv("3", "1", "2", 5000,200, -1, 1, BND, 2);
        final SvVarData dup = createTestSv("4", "1", "1", 1000,6000, -1, 1, DUP, 1);

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


}
