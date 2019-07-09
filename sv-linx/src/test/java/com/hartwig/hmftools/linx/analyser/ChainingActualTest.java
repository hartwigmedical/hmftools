package com.hartwig.hmftools.linx.analyser;

import static com.hartwig.hmftools.common.variant.structural.StructuralVariantType.BND;
import static com.hartwig.hmftools.common.variant.structural.StructuralVariantType.DEL;
import static com.hartwig.hmftools.common.variant.structural.StructuralVariantType.DUP;
import static com.hartwig.hmftools.common.variant.structural.StructuralVariantType.INV;
import static com.hartwig.hmftools.common.variant.structural.StructuralVariantType.SGL;
import static com.hartwig.hmftools.linx.analyser.SvTestHelper.createTestSv;
import static com.hartwig.hmftools.linx.analysis.ClusterAnnotations.ALL_ANNOTATIONS;
import static com.hartwig.hmftools.linx.types.SvArmCluster.ARM_CL_COMPLEX_FOLDBACK;
import static com.hartwig.hmftools.linx.types.SvArmCluster.ARM_CL_DSB;
import static com.hartwig.hmftools.linx.types.SvArmCluster.ARM_CL_FOLDBACK;
import static com.hartwig.hmftools.linx.types.SvArmCluster.ARM_CL_FOLDBACK_DSB;
import static com.hartwig.hmftools.linx.types.SvArmCluster.ARM_CL_ISOLATED_BE;
import static com.hartwig.hmftools.linx.types.SvArmCluster.ARM_CL_TI_ONLY;
import static com.hartwig.hmftools.linx.types.SvArmCluster.getArmClusterData;
import static com.hartwig.hmftools.linx.types.SvCluster.CLUSTER_ANNONTATION_DM;
import static com.hartwig.hmftools.linx.types.SvLinkedPair.ASSEMBLY_MATCH_MATCHED;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertFalse;
import static org.junit.Assert.assertTrue;

import java.util.HashMap;
import java.util.List;
import java.util.Map;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.linx.cn.HomLossEvent;
import com.hartwig.hmftools.linx.types.SvChain;
import com.hartwig.hmftools.linx.types.SvCluster;
import com.hartwig.hmftools.linx.cn.LohEvent;
import com.hartwig.hmftools.linx.types.SvVarData;

import org.junit.Ignore;
import org.junit.Test;

public class ChainingActualTest
{
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
        final SvVarData var2 = createTestSv("78", "3", "0", 25331584, -1, -1, 1, SGL, 12.2, 0, 2.06, 0, 2.12, "");
        final SvVarData var3 = createTestSv("79", "3", "3", 26663922, 26664498, 1, 1, INV, 11.9, 8.0, 3.94, 3.94, 4.26, "");
        final SvVarData var4 = createTestSv("88", "3", "6", 26431918, 26194040, -1, -1, BND, 11.9, 7.3, 3.85, 3.65, 3.77, "");
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

        // check links
        assertEquals(ASSEMBLY_MATCH_MATCHED, var4.getAssemblyMatchType(false));
        assertEquals(ASSEMBLY_MATCH_MATCHED, var5.getAssemblyMatchType(true));

        // should be assembled when assembles from the same breakend are support again
        assertEquals(ASSEMBLY_MATCH_MATCHED, var5.getAssemblyMatchType(false));

        assertEquals(var7.getLinkedPair(false), var8.getLinkedPair(false));
        assertEquals(ASSEMBLY_MATCH_MATCHED, var7.getAssemblyMatchType(false));
        assertEquals(ASSEMBLY_MATCH_MATCHED, var8.getAssemblyMatchType(false));

        assertEquals(var6.getLinkedPair(false), var8.getLinkedPair(true));
        assertEquals(ASSEMBLY_MATCH_MATCHED, var6.getAssemblyMatchType(false));
        assertEquals(ASSEMBLY_MATCH_MATCHED, var8.getAssemblyMatchType(true));

        // check foldbacks
        assertEquals(var1.getFoldbackLink(true), var1.id());
        assertEquals(var3.getFoldbackLink(true), var3.id());
        assertEquals(var4.getFoldbackLink(true), var4.id());
        assertEquals(var6.getFoldbackLink(true), var7.id());
        assertEquals(var7.getFoldbackLink(true), var6.id());

        // check local topology
        final int[] armClusterData = getArmClusterData(cluster);

        assertEquals(8, cluster.getArmClusters().size());
        assertEquals(3, armClusterData[ARM_CL_ISOLATED_BE]);
        assertEquals(2, armClusterData[ARM_CL_TI_ONLY]);
        assertEquals(0, armClusterData[ARM_CL_DSB]);
        assertEquals(3, armClusterData[ARM_CL_FOLDBACK]);
        assertEquals(0, armClusterData[ARM_CL_FOLDBACK_DSB]);
        assertEquals(0, armClusterData[ARM_CL_COMPLEX_FOLDBACK]);

        // check chains
        assertEquals(1, cluster.getChains().size());
        final SvChain chain = cluster.getChains().get(0);
        assertEquals(10, chain.getLinkCount());
    }

    @Test
    public void testActualDoubleMinuteChaining()
    {
        // based on CPCT2100151T chromosome 7 and EGFR AMP

        SvTestHelper tester = new SvTestHelper();
        tester.logVerbose(true);

        tester.Config.RequiredAnnotations = ALL_ANNOTATIONS;

        // tester.Analyser.getChainFinder().setUseNewMethod(false);

        /* SVs

            Id,ChrStart,ChrEnd,PosStart,PosEnd,OrientStart,OrientEnd,Type,CNStart,CNEnd,CNChgStart,CNChgEnd,Ploidy
            78,7,7,54877946,55110179,-1,-1,INV,27.65,27.9,8.02,8.57,7.14
            79,7,7,55145598,55207090,1,1,INV,27.9,20.51,8.2,8.88,7.14
            80,7,7,55200035,55223094,1,-1,DEL,19.71,19.79,8.1,8.17,7.34
            81,7,7,55204831,55293096,-1,1,DUP,20.51,19.79,8.89,8.16,6.86
            82,7,7,55092635,55588244,1,-1,DEL,27.65,19.53,8.31,7.9,7.19
            83,7,7,54832875,55636355,1,-1,DEL,27.85,27.79,8.23,8.26,1.59
            84,7,7,54832269,55636582,-1,1,DUP,27.85,27.79,24.86,24.85,23.21

            SOLUTION: a closed loop:

            chain(0) 0: pair(83 7:55636000:end & 84r 7:55700000:end) Inferred  length(64000)
            chain(0) 1: pair(84r 7:54800000:start & 84r 7:55700000:end) Inferred ADJAC length(900000)
            chain(0) 2: pair(84r 7:54800000:start & 80 7:55200000:start) Inferred ADJAC length(400000)
            chain(0) 3: pair(80 7:55223000:end & 81 7:55293000:end) Inferred ONLY length(70000)
            chain(0) 4: pair(81 7:55204000:start & 79 7:55207000:end) Inferred ADJAC length(3000)
            chain(0) 5: pair(79 7:55145000:start & 78 7:55100000:end) Inferred ADJAC length(45000)
            chain(0) 6: pair(78 7:54877000:start & 82 7:55092000:start) Inferred ONLY length(215000)
            chain(0) 7: pair(82 7:55588000:end & 84 7:55700000:end) Inferred ONLY length(112000)
            chain(0) 8: pair(84 7:54800000:start & 83 7:54832000:start) Inferred ONLY length(32000)
         */

        // merge 5 clusters with varying levels of copy number change (ie replication) from 4 foldbacks
        String chromosome = "7";

        final SvVarData var1 = createTestSv("78", chromosome, chromosome, 54877000, 55100000, -1, -1, INV, 2);
        final SvVarData var2 = createTestSv("79", chromosome, chromosome, 55145000, 55207000, 1, 1, INV, 2);
        final SvVarData var3 = createTestSv("80", chromosome, chromosome, 55200000, 55223000, 1, -1, DEL, 2);
        final SvVarData var5 = createTestSv("81", chromosome, chromosome, 55204000, 55293000, -1, 1, DUP, 2);
        final SvVarData var4 = createTestSv("82", chromosome, chromosome, 55092000, 55588000, 1, -1, DEL, 2);
        final SvVarData var6 = createTestSv("83", chromosome, chromosome, 54832000, 55636000, 1, -1, DEL, 2);
        SvVarData dmDup = createTestSv("84", chromosome, chromosome, 54800000, 55700000, -1, 1, DUP, 10);

        tester.AllVariants.add(dmDup);
        tester.AllVariants.add(var1);
        tester.AllVariants.add(var2);
        tester.AllVariants.add(var3);
        tester.AllVariants.add(var4);
        tester.AllVariants.add(var5);
        tester.AllVariants.add(var6);

        tester.preClusteringInit();

        tester.Analyser.clusterAndAnalyse();

        // check clustering
        assertEquals(tester.Analyser.getClusters().size(), 1);

        final SvCluster cluster = tester.Analyser.getClusters().get(0);
        assertTrue(cluster.getAnnotations().contains(CLUSTER_ANNONTATION_DM));
        assertEquals(1, cluster.getDoubleMinuteSVs().size());

        // check chains
        assertEquals(1, cluster.getChains().size());
        final SvChain chain = cluster.getChains().get(0);
        assertEquals(11, chain.getLinkCount());
        assertTrue(chain.isClosedLoop());
    }


    @Ignore
    @Test
    public void testActualSimpleChaining1()
    {
        // based on CPCT02010325T clusterId 37, where the shortest TI of 76 bases is actually ignored so as to make 2 chains
        SvTestHelper tester = new SvTestHelper();
        // tester.logVerbose(true);
        tester.Analyser.getChainFinder().setMaxPossibleLinks(2);

        final SvVarData var1 = createTestSv("7821420","18","X",23601785,48007145,-1,1,BND,1.92,1.96,0.98,0.95,1.03, "");
        final SvVarData var2 = createTestSv("7821421","X","X",48004021,48123140,-1,-1,INV,0.92,0.99,0.92,0.99,0.96, "");
        final SvVarData var3 = createTestSv("7821422","X","X",48082005,66755692,1,1,INV,1.01,1,1.01,1,0.86, "");
        final SvVarData var4 = createTestSv("7821423","18","X",23577410,66767221,1,-1,BND,1.95,0.97,1.02,0.97,1.01, "");
        final SvVarData var5 = createTestSv("7821424","X","X",47973211,67907761,1,1, INV,1,0.97,1,0.97,0.99, "");
        final SvVarData var6 = createTestSv("7821425","X","X",48007069,67910047,-1,-1,INV,1.96,1,1.04,1,0.99, "");

        tester.AllVariants.add(var1);
        tester.AllVariants.add(var2);
        tester.AllVariants.add(var3);
        tester.AllVariants.add(var4);
        tester.AllVariants.add(var5);
        tester.AllVariants.add(var6);

        List<LohEvent> lohData = tester.CnDataLoader.getLohData();

        lohData.add(new LohEvent( "18", 23601785, 23577410,
                "BND", "BND", 1, 1, 1, 0, 1, 1,
                var1.dbId(), var4.dbId()));

        tester.ClusteringMethods.setSampleCnEventData(lohData, tester.CnDataLoader.getHomLossData());
        tester.preClusteringInit();

        tester.Analyser.clusterAndAnalyse();

        // now check final chain-finding across all sub-clusters
        assertEquals(tester.Analyser.getClusters().size(), 1);
        final SvCluster cluster = tester.Analyser.getClusters().get(0);

        assertEquals(2, cluster.getChains().size());
    }


}
