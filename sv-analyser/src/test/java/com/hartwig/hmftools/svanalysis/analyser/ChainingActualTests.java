package com.hartwig.hmftools.svanalysis.analyser;

import static com.hartwig.hmftools.common.variant.structural.StructuralVariantType.BND;
import static com.hartwig.hmftools.common.variant.structural.StructuralVariantType.DEL;
import static com.hartwig.hmftools.common.variant.structural.StructuralVariantType.DUP;
import static com.hartwig.hmftools.common.variant.structural.StructuralVariantType.INV;
import static com.hartwig.hmftools.common.variant.structural.StructuralVariantType.SGL;
import static com.hartwig.hmftools.svanalysis.analyser.SvTestHelper.createTestSv;
import static com.hartwig.hmftools.svanalysis.types.SvLinkedPair.ASSEMBLY_MATCH_MATCHED;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertTrue;

import java.util.HashMap;
import java.util.List;
import java.util.Map;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.svanalysis.types.SvChain;
import com.hartwig.hmftools.svanalysis.types.SvCluster;
import com.hartwig.hmftools.svanalysis.types.SvLOH;
import com.hartwig.hmftools.svanalysis.types.SvVarData;

import org.junit.Ignore;
import org.junit.Test;

public class ChainingActualTests
{

    @Ignore
    @Test
    public void testActualComplexChaining1()
    {
        // from sampleId CPCT02020258T but not sure if has clustered all SVs correctly
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

        assertEquals(2, cluster.getChains().size());

        assertEquals(14, cluster.getChains().get(0).getLinkCount());
        assertEquals(3, cluster.getChains().get(1).getLinkCount());
    }

    @Test
    public void testActualComplexChaining2()
    {
        // based on COLO829T chromosomes 3 + 6,10,12 and 1

        SvTestHelper tester = new SvTestHelper();
        tester.logVerbose(true);

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

            09:16:20.469 [main] [DEBUG] chain(0) 0: pair(78 3:25331584 SGL-on-known & 119 3:25400602:start) TI inferred len=69018
            09:16:20.469 [main] [DEBUG] chain(0) 1: pair(119 12:72666892:end & 120 12:72667075:end) TI assembly len=183
            09:16:20.470 [main] [DEBUG] chain(0) 2: pair(120 10:60477422:start & 113 10:60477224:end) TI assembly len=198
            09:16:20.470 [main] [DEBUG] chain(0) 3: pair(113 3:25401059:start & 77 3:24566180:end) TI inferred len=834879
            09:16:20.470 [main] [DEBUG] chain(0) 4: pair(77 3:24565108:start & 79 3:26664498:end) TI inferred len=2099390
            09:16:20.470 [main] [DEBUG] chain(0) 5: pair(79 3:26663922:start & 88r 3:26431918:start) TI inferred len=232004
            09:16:20.470 [main] [DEBUG] chain(0) 6: pair(88r 6:26194040:end & 89 6:26194406:end) TI assembly len=366
            09:16:20.470 [main] [DEBUG] chain(0) 7: pair(89 6:26194117:start & 88 6:26194040:end) TI assembly len=77
            09:16:20.470 [main] [DEBUG] chain(0) 8: pair(88 3:26431918:start & 79r 3:26663922:start) TI inferred len=232004
            09:16:20.470 [main] [DEBUG] chain(0) 9: pair(79r 3:26664498:end & 77r 3:24565108:start) TI inferred len=2099390
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
        assertTrue(cluster.hasReplicatedSVs());
        assertEquals(cluster.getSvCount(), 8);
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
        assertEquals(1, cluster.getChains().size());
        final SvChain chain = cluster.getChains().get(0);
        assertEquals(10, chain.getLinkCount());

    }

    @Test
    public void testActualSimpleChaining1()
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
    public void testActualComplexChaining3()
    {
        // based on CPCT02210035, chr 12 with 2 foldbacks, a quasi foldback and a complex DUP (CENTRO - A - B - C - B - D - A - B - C - B - E - 17)
        SvTestHelper tester = new SvTestHelper();
        tester.logVerbose(true);

        final SvVarData var1 = createTestSv("536","12","12",61264340,65981274,-1,1, DUP,3.05,3.03,2.24,2.02,1.95,"");
        final SvVarData var2 = createTestSv("537","12","12",65910212,65983577,1,-1, DEL,2.85,2.94,1.85,1.92,1.72,"");
        final SvVarData var3 = createTestSv("538","12","12",65980999,65990480,-1,1, DUP,3.03,2.94,2.03,1.84,2.06,"");
        final SvVarData var4 = createTestSv("539","12","12",61263455,66021729,1,-1, DEL,2.98,2.96,2.17,2.17,2.17,"");
        final SvVarData var5 = createTestSv("540","12","12",61265162,66158890,1,-1, DEL,3.05,4.98,2.02,2.02,2.02,"");
        final SvVarData var6 = createTestSv("541","12","12",66021681,66159101,-1,1, DUP,3.31,4.98,2.22,1.86,2.06,"");
        final SvVarData var7 = createTestSv("542","12","12",66021726,66160271,1,1, INV,3.31,5.13,2.53,2.17,2.24,"");
        final SvVarData var8 = createTestSv("543","12","12",66159115,66232527,-1,1, DUP,5.13,5.03,2.01,2.05,1.82,"");

        // foldback (event B)
        final SvVarData var9 = createTestSv("545","12","12",68517628,68519581,-1,-1, INV,12.84,20.85,7.85,8.01,7.19,"");

        // complex DUP
        final SvVarData var10 = createTestSv("547","12","12",66303014,69950579,-1,1, DUP,4.97,17.06,2,2.3,1.92,"");
        final SvVarData var11 = createTestSv("548","12","12",65910168,69953088,-1,1, DUP,2.85,14.76,1.79,1.56,1.87,"");

        // inversion (event C)
        final SvVarData var12 = createTestSv("550","12","12",69609346,70421507,1,1, INV,20.85,13.2,3.79,4.25,4.25,"");

        // foldback (event A)
        final SvVarData var13 = createTestSv("551","12","12",70872871,70876753,1,1, INV,8.95,4.78,4.17,3.66,3.74,"");

        final SvVarData var14 = createTestSv("552","4","12",44130417,70876777,1,-1, BND,2.87,2.02,0.89,0.9,0.81,"");
        final SvVarData var15 = createTestSv("603","12","17",66164961,52543282,-1,-1, BND,5.03,4.22,2.07,2.4,1.91,"");
        final SvVarData var16 = createTestSv("604","17","17",52543279,52543541,1,1, INV,3.24,4.22,1.42,1.42,1.45,"");
        final SvVarData var17 = createTestSv("605","17","17",52543194,52543771,-1,-1, INV,3.24,4,1.2,1.2,1.21,"");


        // mark assembled links
        var1.setAssemblyData(false, "asmb13");
        var3.setAssemblyData(true, "asmb13");

        var2.setAssemblyData(true, "asmb211");
        var11.setAssemblyData(true, "asmb211");

        var6.setAssemblyData(true, "asmb67");
        var7.setAssemblyData(true, "asmb67");

        var5.setAssemblyData(false, "asmb56");
        var6.setAssemblyData(false, "asmb56");

        var15.setAssemblyData(false, "asmb1516");
        var16.setAssemblyData(false, "asmb1516");

        var16.setAssemblyData(true, "asmb1617");
        var17.setAssemblyData(true, "asmb1617");

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
        tester.AllVariants.add(var13);
        tester.AllVariants.add(var14);
        tester.AllVariants.add(var15);
        tester.AllVariants.add(var16);
        tester.AllVariants.add(var17);

        tester.preClusteringInit();

        tester.Analyser.clusterAndAnalyse();

        // first check foldbacks are in place
        assertEquals(var9.getFoldbackLink(true), var9.id());
        assertEquals(var13.getFoldbackLink(true), var13.id());

        // now check final chain-finding across all sub-clusters
        assertEquals(tester.Analyser.getClusters().size(), 1);
        final SvCluster cluster = tester.Analyser.getClusters().get(0);

        assertEquals(2, cluster.getChains().size());

        // assertEquals(14, cluster.getChains().get(0).getLinkCount());
        // assertEquals(3, cluster.getChains().get(1).getLinkCount());
    }

    @Test
    @Ignore
    public void testActualComplexChaining4()
    {
        // based on CPCT02080180T with 3 foldbacks
        SvTestHelper tester = new SvTestHelper();
        tester.logVerbose(true);
        tester.Analyser.getChainFinder().setLogVerbose(true);

        SvVarData var1 = createTestSv("900","17","17",63727835,63729894,-1,1,DUP,7.79,5.33,1.04,1.04,1.02,"");
        SvVarData var2 = createTestSv("899","17","17",60677725,63728985,1,1,INV,3.81,5.94,0.69,0.61,0.64,"");
        SvVarData var3 = createTestSv("901","17","17",53291465,63857338,1,1,INV,6.03,5.96,3.01,2.8,2.79,"");
        SvVarData var4 = createTestSv("903","17","17",61869930,66092317,-1,1,DUP,4.19,4.11,1.01,1.03,1.11,"");
        SvVarData var5 = createTestSv("904","17","17",54264135,66294080,-1,-1,INV,4.23,4.13,1.08,1.05,0.88,"");
        SvVarData var6 = createTestSv("905","17","17",63856359,66440216,-1,1,DUP,5.96,6.19,1.67,1.44,1.56,"");
        SvVarData var7 = createTestSv("896","17","17",57423968,60669111,1,1,INV,5.15,5.14,2.05,2.04,2.13,"");
        SvVarData var8 = createTestSv("908","17","17",66436776,66453504,-1,1,DUP,7.47,4.47,1.12,0.24,1.13,"");
        SvVarData var9 = createTestSv("910","17","17",63721643,66759055,-1,1,DUP,6.75,5.3,1.42,0.79,0.93,"");
        SvVarData var10 = createTestSv("912","17","17",65376257,70213096,-1,-1,INV,4.11,4.16,0.95,0.99,1.02,"");
        SvVarData var11 = createTestSv("907","17","17",66435843,66442464,-1,1,DUP,6.35,6.76,2.22,2.29,2.28,"");
        SvVarData var12 = createTestSv("909","17","17",66438835,66755919,1,-1,DEL,7.47,5.3,1.28,1.08,1.28,"");
        SvVarData var13 = createTestSv("911","17","17",53813746,66761458,-1,1,DUP,4.24,4.51,1.12,1.35,0.84,"");
        SvVarData var14 = createTestSv("906","17","17",63728656,66441779,1,-1,DEL,7.79,6.76,1.85,2.01,1.88,"");
        SvVarData var15 = createTestSv("920","17","17",62335817,79024217,1,-1,DEL,4.19,4.13,1.08,1.02,0.86,"");
        SvVarData var16 = createTestSv("913","17","17",60677626,73851024,-1,-1,INV,3.81,3.82,0.71,0.71,0.71,"");
        SvVarData var17 = createTestSv("914","17","17",71119404,73853333,1,-1,DEL,4.16,4.9,1.05,1.08,0.77,"");
        SvVarData var18 = createTestSv("892","17","17",53422106,56798618,1,-1,DEL,5.2,5.15,2.08,2.01,2.03,"");
        SvVarData var19 = createTestSv("893","17","17",53292876,58826153,-1,1,DUP,5.2,5.19,2.18,2.04,2.32,"");
        SvVarData var20 = createTestSv("895","17","17",59899569,60451461,-1,-1,INV,5.14,5.14,1.99,2.05,2.01,"");
        SvVarData var21 = createTestSv("915","17","17",73854033,73858713,1,1,INV,6.12,5.01,1.11,1.11,0.98,"");
        SvVarData var22 = createTestSv("916","17","17",73853710,73858861,-1,1,DUP,6.12,3.9,1.22,0.75,1.04,"");
        SvVarData var23 = createTestSv("917","17","17",73859108,73862099,-1,-1,INV,4.68,6.4,1.53,1.72,1.46,"");
        SvVarData var24 = createTestSv("894","17","17",60447627,60450127,1,1,INV,5.14,4.13,1.01,1.04,1.07,"");
        SvVarData var25 = createTestSv("919","17","17",57571922,75330355,-1,1,DUP,5.19,5.13,2.09,2.02,2.02,"");
        SvVarData var26 = createTestSv("918","17","17",63715443,73862172,-1,1,DUP,5.33,6.4,2.22,1.28,1.33,"");
        SvVarData var27 = createTestSv("889","17","17",53286924,53816913,1,1,INV,5.04,4.24,1.08,1.09,0.95,"");
        SvVarData var28 = createTestSv("888","17","17",53280088,53288151,-1,-1,INV,5.04,6.03,1.92,2.07,2.04,"");
        SvVarData var29 = createTestSv("891","17","17",51738130,54666910,1,1,INV,4.15,4.17,1.03,1.03,1.03,"");

        var2.setAssemblyData(true, "asmb_2_16");
        var16.setAssemblyData(true, "asmb_2_16");

        var23.setAssemblyData(false, "asmb_23_26");
        var26.setAssemblyData(false, "asmb_23_26");

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
        tester.AllVariants.add(var13);
        tester.AllVariants.add(var14);
        tester.AllVariants.add(var15);
        tester.AllVariants.add(var16);
        tester.AllVariants.add(var17);
        tester.AllVariants.add(var18);
        tester.AllVariants.add(var19);
        tester.AllVariants.add(var20);
        tester.AllVariants.add(var21);
        tester.AllVariants.add(var22);
        tester.AllVariants.add(var23);
        tester.AllVariants.add(var24);
        tester.AllVariants.add(var25);
        tester.AllVariants.add(var26);
        tester.AllVariants.add(var27);
        tester.AllVariants.add(var28);
        tester.AllVariants.add(var29);

        tester.preClusteringInit();

        tester.Analyser.clusterAndAnalyse();

        // first check foldbacks are in place
        assertEquals(var21.getFoldbackLink(true), var21.id());
        assertEquals(var23.getFoldbackLink(true), var23.id());
        assertEquals(var24.getFoldbackLink(true), var24.id());

        // now check final chain-finding across all sub-clusters
        assertEquals(tester.Analyser.getClusters().size(), 1);
        final SvCluster cluster = tester.Analyser.getClusters().get(0);

        // assertEquals(cluster.getChains().size(), 6);

    }


}
