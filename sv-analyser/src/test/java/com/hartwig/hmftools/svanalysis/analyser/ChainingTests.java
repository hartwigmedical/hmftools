package com.hartwig.hmftools.svanalysis.analyser;

import static com.hartwig.hmftools.common.variant.structural.StructuralVariantType.BND;
import static com.hartwig.hmftools.common.variant.structural.StructuralVariantType.DEL;
import static com.hartwig.hmftools.common.variant.structural.StructuralVariantType.DUP;
import static com.hartwig.hmftools.common.variant.structural.StructuralVariantType.INV;
import static com.hartwig.hmftools.common.variant.structural.StructuralVariantType.SGL;
import static com.hartwig.hmftools.svanalysis.analyser.SvTestHelper.createBnd;
import static com.hartwig.hmftools.svanalysis.analyser.SvTestHelper.createDel;
import static com.hartwig.hmftools.svanalysis.analyser.SvTestHelper.createDup;
import static com.hartwig.hmftools.svanalysis.analyser.SvTestHelper.createInv;
import static com.hartwig.hmftools.svanalysis.analyser.SvTestHelper.createTestSv;
import static com.hartwig.hmftools.svanalysis.types.SvChain.CHAIN_ASSEMBLY_LINK_COUNT;
import static com.hartwig.hmftools.svanalysis.types.SvChain.CHAIN_LINK_COUNT;
import static com.hartwig.hmftools.svanalysis.types.SvLinkedPair.ASSEMBLY_MATCH_MATCHED;
import static com.hartwig.hmftools.svanalysis.types.SvLinkedPair.LINK_TYPE_DB;
import static com.hartwig.hmftools.svanalysis.types.SvLinkedPair.LINK_TYPE_TI;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertFalse;
import static org.junit.Assert.assertTrue;

import java.util.HashMap;
import java.util.List;
import java.util.Map;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.svanalysis.types.SvChain;
import com.hartwig.hmftools.svanalysis.types.SvCluster;
import com.hartwig.hmftools.svanalysis.types.SvLOH;
import com.hartwig.hmftools.svanalysis.types.SvLinkedPair;
import com.hartwig.hmftools.svanalysis.types.SvVarData;

import org.junit.Ignore;
import org.junit.Test;

public class ChainingTests
{

    @Test
    public void testBasicChainBuilding()
    {
        // create a chain out of simple DELs and test the various chaining features
        final SvVarData var1 = createDel("1", "1", 1100, 1200);
        final SvVarData var2 = createDel("2", "1", 1300, 1400);
        SvLinkedPair lp1 = new SvLinkedPair(var1, var2, LINK_TYPE_TI, false, true);
        lp1.setIsInferred(false);

        // test adding linked pairs of various orientations to the start and end of a chain
        SvChain chain = new SvChain(0);

        chain.addLink(lp1, true);

        assertTrue(chain.firstLinkOpenOnStart());
        assertTrue(!chain.lastLinkOpenOnStart());

        final SvVarData var3 = createDel("3", "1", 1500, 1600);
        SvLinkedPair lp2 = new SvLinkedPair(var3, var2, LINK_TYPE_TI, true, false);
        lp2.setIsInferred(false);

        assertFalse(chain.canAddLinkedPairToStart(lp2));
        assertTrue(chain.canAddLinkedPairToEnd(lp2));

        chain.addLink(lp2, false);

        assertEquals(chain.getFirstSV(), var1);
        assertEquals(chain.getLastSV(), var3);
        assertEquals(chain.getSvList().get(1), var2);

        final SvVarData var4 = createDel("4", "1", 900, 1000);
        SvLinkedPair lp3 = new SvLinkedPair(var1, var4, LINK_TYPE_TI, true, false);

        assertTrue(chain.canAddLinkedPairToStart(lp3));
        assertFalse(chain.canAddLinkedPairToEnd(lp3));
        chain.addLink(lp3, true);

        assertEquals(chain.getFirstSV(), var4);
        assertEquals(chain.getLastSV(), var3);
        assertEquals(chain.getSvList().get(1), var1);
        assertEquals(chain.getSvList().get(2), var2);

        // test a potentially closing link
        final SvVarData var5 = createDup("5", "1", 800, 1700);
        SvLinkedPair lp4 = new SvLinkedPair(var5, var4, LINK_TYPE_TI, true, true);

        assertTrue(chain.canAddLinkedPairToStart(lp4));
        chain.addLink(lp4, true);

        assertEquals(chain.getAssemblyLinkCount(), 2);

        SvLinkedPair lp5 = new SvLinkedPair(var5, var3, LINK_TYPE_TI, false, false);

        assertTrue(chain.canAddLinkedPairToEnd(lp5));
        assertTrue(chain.linkWouldCloseChain(lp5));

        // tests paths through the chain from various points
        int[] chainData = chain.breakendsAreChained(var4, false, var3, true);
        assertEquals(chainData[CHAIN_LINK_COUNT], 3);
        assertEquals(chainData[CHAIN_ASSEMBLY_LINK_COUNT], 2);

        // check works in the other direction
        chainData = chain.breakendsAreChained(var3, true, var4, false);
        assertEquals(chainData[CHAIN_LINK_COUNT], 3);
        assertEquals(chainData[CHAIN_ASSEMBLY_LINK_COUNT], 2);

        // check a single link
        chainData = chain.breakendsAreChained(var1, false, var2, true);
        assertEquals(chainData[CHAIN_LINK_COUNT], 1);
        assertEquals(chainData[CHAIN_ASSEMBLY_LINK_COUNT], 1);

        // check breakends facing the wrong way
        chainData = chain.breakendsAreChained(var1, false, var2, false);
        assertEquals(chainData[CHAIN_LINK_COUNT], 0);

        chainData = chain.breakendsAreChained(var1, true, var2, false);
        assertEquals(chainData[CHAIN_LINK_COUNT], 0);

        chainData = chain.breakendsAreChained(var1, true, var2, true);
        assertEquals(chainData[CHAIN_LINK_COUNT], 0);

        // check no link
        chainData = chain.breakendsAreChained(var5, false, var1, true);
        assertEquals(chainData[CHAIN_LINK_COUNT], 0);
    }

    @Test
    public void testFullyAssembledChain()
    {
        SvTestHelper tester = new SvTestHelper();
        tester.logVerbose(true);

        tester.Analyser.getChainFinder().setUseNewMethod(true);

        final SvVarData var1 = createDel("0", "1", 100,200);
        final SvVarData var2 = createDel("1", "1", 300,400);
        final SvVarData var3 = createDel("2", "1", 500,600);
        final SvVarData var4 = createDel("3", "1", 700,800);

        var1.setAssemblyData(false, "asmb12");
        var2.setAssemblyData(true, "asmb12");
        var2.setAssemblyData(false, "asmb23");
        var3.setAssemblyData(true, "asmb23");
        var3.setAssemblyData(false, "asmb34");
        var4.setAssemblyData(true, "asmb34");

        // add them out of order which will require partial chain reconciliation
        tester.AllVariants.add(var1);
        tester.AllVariants.add(var2);
        tester.AllVariants.add(var3);
        tester.AllVariants.add(var4);

        tester.preClusteringInit();
        tester.Analyser.clusterAndAnalyse();

        assertEquals(1, tester.Analyser.getClusters().size());
        final SvCluster cluster = tester.Analyser.getClusters().get(0);

        assertEquals(1, cluster.getChains().size());

        final SvChain chain = cluster.getChains().get(0);

        assertEquals(3, chain.getLinkCount());
    }

    @Test
    public void testPartiallyAssembledChain()
    {
        SvTestHelper tester = new SvTestHelper();
        tester.logVerbose(true);

        tester.Analyser.getChainFinder().setUseNewMethod(true);

        final SvVarData var0 = createDel("0", "1", 100,200);
        final SvVarData var1 = createDel("1", "1", 300,400);
        final SvVarData var2 = createDel("2", "1", 500,600);
        final SvVarData var3 = createDel("3", "1", 700,800);

        var1.setAssemblyData(false, "asmb23");
        var2.setAssemblyData(true, "asmb23");

        // add them out of order which will require partial chain reconciliation
        tester.AllVariants.add(var0);
        tester.AllVariants.add(var1);
        tester.AllVariants.add(var2);
        tester.AllVariants.add(var3);

        tester.preClusteringInit();
        tester.Analyser.clusterAndAnalyse();

        assertEquals(1, tester.Analyser.getClusters().size());
        final SvCluster cluster = tester.Analyser.getClusters().get(0);

        assertEquals(1, cluster.getChains().size());

        final SvChain chain = cluster.getChains().get(0);

        assertEquals(3, chain.getLinkCount());
    }

    @Test
    public void testComplexChaining1()
    {
        SvTestHelper tester = new SvTestHelper();
        tester.logVerbose(true);

        final SvVarData var1 = createTestSv("717", "11", "11", 67149357,67150121, -1, -1, INV, 10.4, 18, 7.51, 7.51, 6, "");
        final SvVarData var2 = createTestSv("719", "11", "11", 67465724,68587766, -1, 1, DUP, 20.7, 19.7, 2.75, 1.54, 2, "");
        final SvVarData var3 = createTestSv("731", "11", "11", 68574864,70107010, 1, -1, DEL, 19.7, 21.0, 1.6, 2.05, 1.9, "");
        final SvVarData var4 = createTestSv("728", "11", "11", 68587010,69911892, -1, -1, INV, 19.7, 15.3, 1.67, 1.66, 1.54, "");
        final SvVarData var5 = createTestSv("723", "11", "11", 68797542,68897549, 1, 1, INV, 18.7, 15.6, 3.04, 3.96, 3.33, "");
        final SvVarData var6 = createTestSv("726", "11", "11", 69720626,69722514, -1, -1, INV, 9.8, 13.8, 4.13, 3.94, 4.2, "");
        final SvVarData var7 = createTestSv("727", "11", "11", 69722658,69724290, 1, -1, DEL, 13.8,13.7, 3.98, 3.88, 4.36, "");
        final SvVarData var8 = createTestSv("732", "11", "11", 70106999,70107331, 1, 1, INV, 20.8, 21, 1.85, 1.85, 1.82, "");
        final SvVarData var9 = createTestSv("733", "11", "11", 70170579,70173624, 1, -1, DEL, 19.2, 19.4, 0.67, 0.86, 0.94, "");
        final SvVarData var10 = createTestSv("734", "11", "11", 70173632,70174021, 1, 1, INV, 19.4, 18.3, 1.07, 0.71, 0.78, "");
        final SvVarData var11 = createTestSv("720", "1", "11", 100499653,68672544, 1, -1, BND, 2.05, 19.6, 1, 1.4, 1.2, "");
        final SvVarData var12 = createTestSv("462", "1", "1", 113342526,113343275, -1, -1, INV, 3, 3.9, 0.93, 0.93, 1.11, "");

        // mark assembled links
        var3.setAssemblyData(false, "asmb1");
        var8.setAssemblyData(false, "asmb1");

        var9.setAssemblyData(false, "asmb2");
        var10.setAssemblyData(false, "asmb2");

        var6.setAssemblyData(false, "asmb3");
        var7.setAssemblyData(true, "asmb3");

        tester.AllVariants.add(var1);
        tester.AllVariants.add(var2);
        tester.AllVariants.add(var3);
        tester.AllVariants.add(var4);
        tester.AllVariants.add(var5);
        tester.AllVariants.add(var6);
        tester.AllVariants.add(var7);
        tester.AllVariants.add(var8);
        tester.AllVariants.add(var9);
        tester.AllVariants.add(var10);
        tester.AllVariants.add(var11);
        tester.AllVariants.add(var12);

        Map<String, List<SvLOH>> lohDataMap = new HashMap();
        List<SvLOH> lohData = Lists.newArrayList();

        lohData.add(new SvLOH(tester.SampleId, "1", 1, 2, 100499653, 113342526,
                "BND", "INV", 1, 1, 1, 0, 1, 1,
                "720", "462", false, true));

        lohDataMap.put(tester.SampleId, lohData);

        tester.ClusteringMethods.setSampleLohData(lohDataMap);

        tester.preClusteringInit();

        tester.Analyser.clusterAndAnalyse();

        // first check foldbacks are in place
        assertEquals(var1.getFoldbackLink(true), var1.id());
        assertEquals(var5.getFoldbackLink(true), var5.id());
        assertEquals(var6.getFoldbackLink(true), var6.id());
        assertEquals(var9.getFoldbackLink(true), var10.id());
        assertEquals(var10.getFoldbackLink(true), var9.id());
        assertEquals(var12.getFoldbackLink(true), var12.id());

        // now check final chain-finding across all sub-clusters
        assertEquals(tester.Analyser.getClusters().size(), 1);
        final SvCluster cluster = tester.Analyser.getClusters().get(0);

        assertEquals(cluster.getChains().size(), 2);

        assertEquals(14, cluster.getChains().get(0).getLinkCount());
        assertEquals(3, cluster.getChains().get(1).getLinkCount());
    }

    @Test
    public void testComplexChaining2()
    {
        // based on COLO829T chromosomes 3 + 6,10,12 and 1

        SvTestHelper tester = new SvTestHelper();
        tester.logVerbose(true);
        tester.Analyser.getChainFinder().setUseNewMethod(true);

            /*
            Id	Type	Ploidy	ChrStart	PosStart	OS	AS	CNStart	CNChgS	ChrEnd	PosEnd	    OE	AE	CNEnd	CNChgEnd
            77	INV	    3.96	3	        24565108	-1	P	6.07	4.07	3	    24566180	-1	P	10.14	4.07
            78	SGL	    2.62	3	        25331584	-1	P	12.17	2.03	0	    -1	        0	P	0	    0
            79	INV	    5.26	3	        26663922	1	P	11.92	3.94	3	    26664498	1	P	7.98	3.94
            88	BND	    3.27	3	        26431918	-1	P	11.92	3.83	6	    26194040	-1	P	7.31	3.41
            89	INV	    1.5	    6	        26194117	1	P	7.31	2.04	6	    26194406	1	P	5.27	1.43
            113	BND	    1.77	3	        25401059	1	P	9.94	1.86	10	    60477224	-1	Q	4.06	2.06
            119	BND	    2.19	3	        25400602	1	P	12.17	2.22	12	    72666892	-1	Q	5.22	2.19
            120	BND	    2.26	10	        60477422	1	Q	4.06	2.04	12	    72667075	1	Q	5.22	2.22

            SOLUTION:

            21:30:03.844 [main] [DEBUG] cluster(5) adding complete chain(0) with 10 linked pairs:
            21:30:03.844 [main] [DEBUG] chain(0) 0: pair(78 3:25331584 SGL-on-known & 119 3:25400602:start) TI inferred len=69018
            21:30:03.844 [main] [DEBUG] chain(0) 1: pair(119 12:72666892:end & 120 12:72667075:end) TI assembly len=183
            21:30:03.844 [main] [DEBUG] chain(0) 2: pair(120 10:60477422:start & 113 10:60477224:end) TI assembly len=198
            21:30:03.844 [main] [DEBUG] chain(0) 3: pair(113 3:25401059:start & 77 3:24566180:end) TI inferred len=834879
            21:30:03.845 [main] [DEBUG] chain(0) 4: pair(77 3:24565108:start & 79 3:26663922:start) TI inferred len=2098814
            21:30:03.845 [main] [DEBUG] chain(0) 5: pair(79 3:26664498:end & 88 3:26431918:start) TI inferred len=232580
            21:30:03.845 [main] [DEBUG] chain(0) 6: pair(88 6:26194040:end & 89 6:26194117:start) TI assembly len=77
            21:30:03.845 [main] [DEBUG] chain(0) 7: pair(89 6:26194406:end & 88r 6:26194040:end) TI assembly len=366
            21:30:03.845 [main] [DEBUG] chain(0) 8: pair(88r 3:26431918:start & 79r 3:26663922:start) TI inferred len=232004
            21:30:03.846 [main] [DEBUG] chain(0) 9: pair(79r 3:26664498:end & 77r 3:24565108:start) TI inferred len=2099390
             */

        // merge 5 clusters with varying levels of copy number change (ie replication) from 4 foldbacks
        final SvVarData var1 = createTestSv("77", "3", "3", 24565108, 24566180, -1, -1, INV, 6.1, 10.1, 4.07, 4.07, 3.96, "");
        final SvVarData var2 = createTestSv("78", "3", "0", 25331584, -1, -1, 1, SGL, 12.2, 0, 2.06, 0, 2.62, "");
        final SvVarData var3 = createTestSv("79", "3", "3", 26663922, 26664498, 1, 1, INV, 11.9, 8.0, 3.94, 3.94, 5.26, "");
        final SvVarData var4 = createTestSv("88", "3", "6", 26431918, 26194040, -1, -1, BND, 11.9, 7.3, 3.85, 3.35, 3.27, "");
        final SvVarData var5 = createTestSv("89", "6", "6", 26194117, 26194406, 1, 1, INV, 7.3, 5.3, 2.06, 1.43, 1.5, "");
        final SvVarData var6 = createTestSv("113", "3", "10", 25401059, 60477224, 1, -1, BND, 9.9, 4.1, 1.92, 2.02, 1.77, "");
        final SvVarData var7 = createTestSv("119", "3", "12", 25400602, 72666892, 1, -1, BND, 12.2, 5.2, 2.22, 2.18, 2.19, "");
        final SvVarData var8 = createTestSv("120", "10", "12", 60477422, 72667075, 1, 1, BND, 4.1, 5.2, 2.01, 2.16, 2.26, "");

        // mark assembled links
        var4.setAssemblyData(false, "asmb1a;asmb1b");
        var5.setAssemblyData(true, "asmb1a");
        var5.setAssemblyData(false, "asmb1b");

        var7.setAssemblyData(false, "asmb2");
        var8.setAssemblyData(false, "asmb2");

        var6.setAssemblyData(false, "asmb3");
        var8.setAssemblyData(true, "asmb3");

        // cluster
        tester.AllVariants.add(var1);
        tester.AllVariants.add(var2);
        tester.AllVariants.add(var3);
        tester.AllVariants.add(var4);
        tester.AllVariants.add(var5);
        tester.AllVariants.add(var6);
        tester.AllVariants.add(var7);
        tester.AllVariants.add(var8);

        tester.preClusteringInit();

        tester.Analyser.clusterAndAnalyse();

        // check clustering
        assertEquals(tester.Analyser.getClusters().size(), 1);

        final SvCluster cluster = tester.Analyser.getClusters().get(0);
        assertEquals(cluster.getUniqueSvCount(), 8);
        assertEquals(cluster.getSubClusters().size(), 5);

        // check links
        assertEquals(var4.getAssemblyMatchType(false), ASSEMBLY_MATCH_MATCHED);
        assertEquals(var5.getAssemblyMatchType(true), ASSEMBLY_MATCH_MATCHED);
        assertEquals(var5.getAssemblyMatchType(false), ASSEMBLY_MATCH_MATCHED);

        assertEquals(var7.getLinkedPair(false), var8.getLinkedPair(false));
        assertEquals(var7.getAssemblyMatchType(false), ASSEMBLY_MATCH_MATCHED);
        assertEquals(var8.getAssemblyMatchType(false), ASSEMBLY_MATCH_MATCHED);

        assertEquals(var6.getLinkedPair(false), var8.getLinkedPair(true));
        assertEquals(var6.getAssemblyMatchType(false), ASSEMBLY_MATCH_MATCHED);
        assertEquals(var8.getAssemblyMatchType(true), ASSEMBLY_MATCH_MATCHED);

        // check foldbacks
        assertEquals(var1.getFoldbackLink(true), var1.id());
        assertEquals(var3.getFoldbackLink(true), var3.id());
        assertEquals(var6.getFoldbackLink(true), var7.id());
        assertEquals(var7.getFoldbackLink(true), var6.id());

        // check chains
        assertEquals(cluster.getChains().size(), 1);
        final SvChain chain = cluster.getChains().get(0);
        assertEquals(chain.getLinkCount(), 10);

    }

    @Test
    public void testSimpleChaining1()
    {
        // based on CPCT02010325T clusterId 37, where the shortest TI of 76 bases is actually ignored so as to make 2 chains
        SvTestHelper tester = new SvTestHelper();
        tester.logVerbose(true);

        final SvVarData var1 = createTestSv("7821420","18","X",23601785,48007145,-1,1,BND,1.92,1.96,0.98,0.95,1.03, "");
        final SvVarData var2 = createTestSv("7821421","X","X",48004021,48123140,-1,-1,INV,0.92,0.99,0.92,0.99,0.96, "");
        final SvVarData var3 = createTestSv("7821422","X","X",48082005,66755692,1,1,INV,1.01,1,1.01,1,0.86, "");
        final SvVarData var4 = createTestSv("7821423","18","X",23577410,66767221,1,-1,BND,1.95,0.97,1.02,0.97,1.01, "");
        final SvVarData var5 = createTestSv("7821424","X","X",47973211,67907761,1,1,INV,1,0.97,1,0.97,0.99, "");
        final SvVarData var6 = createTestSv("7821425","X","X",48007069,67910047,-1,-1,INV,1.96,1,1.04,1,0.99, "");

        tester.AllVariants.add(var1);
        tester.AllVariants.add(var2);
        tester.AllVariants.add(var3);
        tester.AllVariants.add(var4);
        tester.AllVariants.add(var5);
        tester.AllVariants.add(var6);

        Map<String, List<SvLOH>> lohDataMap = new HashMap();
        List<SvLOH> lohData = Lists.newArrayList();

        lohData.add(new SvLOH(tester.SampleId, "18", 1, 2, 23601785, 23577410,
                "BND", "BND", 1, 1, 1, 0, 1, 1,
                var1.id(), var4.id(), false, true));

        lohDataMap.put(tester.SampleId, lohData);

        tester.ClusteringMethods.setSampleLohData(lohDataMap);

        tester.preClusteringInit();

        tester.Analyser.clusterAndAnalyse();

        // now check final chain-finding across all sub-clusters
        assertEquals(tester.Analyser.getClusters().size(), 1);
        final SvCluster cluster = tester.Analyser.getClusters().get(0);

        assertEquals(2, cluster.getChains().size());

        // assertEquals(14, cluster.getChains().get(0).getLinkCount());
        // assertEquals(3, cluster.getChains().get(1).getLinkCount());

    }

    @Test
    @Ignore
    public void testComplexChaining3()
    {
        SvTestHelper tester = new SvTestHelper();
        tester.logVerbose(true);

        final SvVarData var1 = createTestSv("8349725","1","1",9861138,9861324,1,1,INV,6.06,5.74,1.59,1.59,1.58, "");
        final SvVarData var2 = createTestSv("8349723","1","1",9004324,9005765,1,1,INV,7.96,6.06,1.9,1.9,1.53, "");
        final SvVarData var3 = createTestSv("8349854","1","4",24653435,1801105,1,1,BND,36.28,9.6,4.92,4.51,4.61, "");
        final SvVarData var4 = createTestSv("8349921","1","5",9861144,165461529,-1,1,BND,5.74,4.98,1.27,2.23,1.75, "");
        final SvVarData var5 = createTestSv("8349732","1","1",25580398,25580626,1,1,INV,31.36,28.25,3.11,3.11,1.42, "");
        final SvVarData var6 = createTestSv("8349735","1","1",26014233,26017546,1,1,INV,22.14,16.21,5.93,5.41,5.99, "");
        final SvVarData var7 = createTestSv("8349736","1","1",26064778,26066088,1,1,INV,10.81,8.36,2.44,3.34,2.85, "");
        final SvVarData var8 = createTestSv("8350322","1","15",23973179,95331774,-1,-1,BND,16.79,12.22,14.7,8.24,9.27, "");
        final SvVarData var9 = createTestSv("8350323","1","15",23974263,95331935,-1,1,BND,30.7,12.22,13.91,8.16,8.73, "");
        final SvVarData var10 = createTestSv("8350409","1","19",23966804,12257497,1,1,BND,4.17,6.08,2.08,1.99,2.01, "");
        final SvVarData var11 = createTestSv("8350412","1","19",8718959,34828900,-1,-1,BND,8.29,8.07,2.7,1.89,2.02,"");

        // mark assembled links
        var1.setAssemblyData(false, "asmb1");
        var4.setAssemblyData(true, "asmb1");

        tester.AllVariants.add(var1);
        tester.AllVariants.add(var2);
        tester.AllVariants.add(var3);
        tester.AllVariants.add(var4);
        tester.AllVariants.add(var5);
        tester.AllVariants.add(var6);
        tester.AllVariants.add(var7);
        tester.AllVariants.add(var8);
        tester.AllVariants.add(var9);
        tester.AllVariants.add(var10);
        tester.AllVariants.add(var11);

        Map<String, List<SvLOH>> lohDataMap = new HashMap();
        List<SvLOH> lohData = Lists.newArrayList();

        lohData.add(new SvLOH(tester.SampleId, "1", 1, 2, 23973179, 23966804,
                "BND", "BND", 1, 1, 1, 0, 1, 1,
                var8.id(), var10.id(), false, true));

        lohDataMap.put(tester.SampleId, lohData);

        tester.ClusteringMethods.setSampleLohData(lohDataMap);

        tester.preClusteringInit();

        tester.Analyser.clusterAndAnalyse();

        // first check foldbacks are in place
        assertEquals(var2.getFoldbackLink(true), var2.id());
        assertEquals(var5.getFoldbackLink(true), var5.id());
        assertEquals(var6.getFoldbackLink(true), var6.id());
        assertEquals(var5.getFoldbackLink(true), var7.id());
        assertEquals(var8.getFoldbackLink(true), var9.id());
        assertEquals(var9.getFoldbackLink(true), var8.id());

        // now check final chain-finding across all sub-clusters
        assertEquals(tester.Analyser.getClusters().size(), 1);
        final SvCluster cluster = tester.Analyser.getClusters().get(0);

        assertEquals(cluster.getChains().size(), 2);

        assertEquals(14, cluster.getChains().get(0).getLinkCount());
        assertEquals(3, cluster.getChains().get(1).getLinkCount());
    }

}
