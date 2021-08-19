package com.hartwig.hmftools.linx.chaining;

import static com.hartwig.hmftools.linx.analysis.ClusterAnnotations.ALL_ANNOTATIONS;
import static com.hartwig.hmftools.linx.analysis.ClusteringReason.FOLDBACKS;
import static com.hartwig.hmftools.linx.chaining.ChainUtils.identicalChain;
import static com.hartwig.hmftools.linx.types.ArmCluster.ARM_CL_COMPLEX_FOLDBACK;
import static com.hartwig.hmftools.linx.types.ArmCluster.ARM_CL_DSB;
import static com.hartwig.hmftools.linx.types.ArmCluster.ARM_CL_FOLDBACK;
import static com.hartwig.hmftools.linx.types.ArmCluster.ARM_CL_FOLDBACK_DSB;
import static com.hartwig.hmftools.linx.types.ArmCluster.ARM_CL_ISOLATED_BE;
import static com.hartwig.hmftools.linx.types.ArmCluster.ARM_CL_TI_ONLY;
import static com.hartwig.hmftools.linx.types.ArmCluster.getArmClusterData;
import static com.hartwig.hmftools.linx.types.SvCluster.CLUSTER_ANNOT_DM;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertTrue;

import java.util.List;

import com.hartwig.hmftools.linx.types.SvBreakend;
import com.hartwig.hmftools.linx.types.SvCluster;
import com.hartwig.hmftools.linx.cn.LohEvent;
import com.hartwig.hmftools.linx.types.SvVarData;

import org.junit.Test;

import com.hartwig.hmftools.linx.utils.LinxTester;
import com.hartwig.hmftools.linx.utils.SampleDataLoader;

// tests modelled on examples from actual samples
public class ChainingActualTest
{
    @Test
    public void testActualComplexChaining()
    {
        // based on COLO829T chromosomes 3 + 6,10,12 and 1

        LinxTester tester = new LinxTester();

        // Configurator.setRootLevel(Level.DEBUG);

        final List<SvVarData> svList = SampleDataLoader.loadSampleTestData("COLO829T");

        final SvVarData var1 = svList.get(0);
        final SvVarData var2 = svList.get(1);
        final SvVarData var3 = svList.get(2);
        final SvVarData var4 = svList.get(3);
        final SvVarData var5 = svList.get(4);
        final SvVarData var6 = svList.get(5);
        final SvVarData var7 = svList.get(6);
        final SvVarData var8 = svList.get(7);

        tester.AllVariants.addAll(svList);

        tester.preClusteringInit();

        tester.Analyser.clusterAndAnalyse();

        // check clustering
        assertEquals(tester.Analyser.getClusters().size(), 1);

        final SvCluster cluster = tester.Analyser.getClusters().get(0);
        assertTrue(cluster.requiresReplication());
        assertEquals(8, cluster.getSvCount());

        // check links
        assertTrue(var4.hasAssemblyLink(false));
        assertTrue(var5.hasAssemblyLink(true));

        // should be assembled when assembles from the same breakend are support again
        assertTrue(var5.hasAssemblyLink(false));

        assertEquals(var7.getLinkedPair(false), var8.getLinkedPair(false));
        assertTrue(var7.hasAssemblyLink(false));
        assertTrue(var8.hasAssemblyLink(false));

        assertEquals(var6.getLinkedPair(false), var8.getLinkedPair(true));
        assertTrue(var6.hasAssemblyLink(false));
        assertTrue(var8.hasAssemblyLink(true));

        // check foldbacks
        assertEquals(var1.getFoldbackId(true), var1.id());
        assertEquals(var3.getFoldbackId(true), var3.id());
        assertEquals(var4.getFoldbackId(true), var4.id());
        assertEquals(var6.getFoldbackId(true), var7.id());
        assertEquals(var7.getFoldbackId(true), var6.id());

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
        assertEquals(8, chain.getSvCount());

        final SvBreakend chainStart = chain.getOpenBreakend(true);
        final SvBreakend chainEnd = chain.getOpenBreakend(false);

        assertTrue((chainStart == null && chainEnd.getSV() == var1) || (chainEnd == null && chainStart.getSV() == var1));

        /* strip out and store log [0-9][0-9]:[0-9][0-9]:[0-9][0-9].[0-9][0-9][0-9] \[main\] \[[A-Z].*\]

         cluster(0) adding assembly linked TI pair(119 12:72666892:end & 120 12:72667075:end) length(183)
         cluster(0) adding assembly linked TI pair(113 10:60477224:end & 120 10:60477422:start) length(198)
         index(0) method(ASSEMBLY) adding linked pair(119 12:72666892:end & 120 12:72667075:end ploidy=1.0) to new chain(0) ploidy(1.0 unc=0.5)
         end breakend exhausted: id(119) ploidy(1.0-1.0-1.0) counts(s=0.0 e=1.0)
         end breakend exhausted: id(120) ploidy(1.0-1.0-1.0) counts(s=0.0 e=1.0)
         pair(113 10:60477224:end & 120 10:60477422:start) links second breakend to chain(0) as only exhausted connection
         index(1) method(ASSEMBLY) adding linked pair(113 10:60477224:end & 120 10:60477422:start ploidy=1.0) to existing chain(0) ploidy(1.0 unc=0.4)
         end breakend exhausted: id(113) ploidy(1.0-1.0-1.0) counts(s=0.0 e=1.0)
         start breakend exhausted: id(120) ploidy(1.0-1.0-1.0) counts(s=1.0 e=1.0)
         SV(id(120) ploidy(1.0-1.0-1.0) counts(s=1.0 e=1.0)) both breakends exhausted
         cluster(0) chaining finished: chains(1 links=2) SVs(3) unlinked SVs(0) breakends(2)
         cluster(0) added chain(0) ploidy(1.0) with 2 linked pairs:
         chain(0): 3_P_T - s_119_e - e_120_s - e_113_s - 3_P_T
         chain(0) 0: pair(119 12:72666892:end & 120 12:72667075:end) ASSEMBLY length(183) index(0)
         chain(0) 1: pair(120 10:60477422:start & 113 10:60477224:end) ASSEMBLY length(198) index(1)
         cluster(3) adding assembly linked TI pair(88 6:26194040:end & 89 6:26194117:start) length(77)
         cluster(3) adding assembly linked TI pair(88 6:26194040:end & 89 6:26194406:end) length(366)
         cluster(3) SV(id(88) pos(3:-1:26431918 -> 6:-1:26194040)) ploidy multiple(2, ploidy(4 vs min=2.0)
         assembly multi-sgl-conn pair(88 6:26194040:end & 89 6:26194117:start) ploidy(1.5): first(ploidy=3.5 links=2) second(ploidy=1.5 links=1)
         index(0) method(ASSEMBLY) adding linked pair(88 6:26194040:end & 89 6:26194117:start ploidy=1.5) to new chain(0) ploidy(1.5 unc=0.6)
         assembly multi-sgl-conn pair(88 6:26194040:end & 89 6:26194406:end) ploidy(1.8): first(ploidy=2.0 links=2) second(ploidy=1.5 links=1)
         pair(88 6:26194040:end & 89 6:26194406:end) links second breakend to chain(0) as only exhausted connection
         index(1) method(ASSEMBLY) adding linked pair(88 6:26194040:end & 89 6:26194406:end ploidy=1.8) to existing chain(0) ploidy(1.7 unc=0.5)
         created 1 partial chains from 2 assembly links
         cluster(3) chaining finished: chains(1 unique=1 links=2) SVs(2) unlinked SVs(0 ploidy=3.5) breakends(1 ploidy=3.7)
         cluster(3) added chain(0) ploidy(1.7) with 2 linked pairs:
         chain(0): 3_P_C - s_88_e - s_89_e - e_88_s - 3_P_C
         chain(0) 0: pair(88 6:26194040:end & 89 6:26194117:start) ASSEMBLY length(77) index(0)
         chain(0) 1: pair(89 6:26194406:end & 88 6:26194040:end) ASSEMBLY length(366) index(1)
         cluster(0) complex SVs(3) desc(BND=3 res=NONE) arms(3) consis(2) chains(1 perc=1.00) replic(false)
         cluster(3) complex SVs(2) desc(BND=1_INV=1 res=NONE) arms(2) consis(0) chains(1 perc=1.00) replic(true) inv=1
         cluster(1) foldback inversion SV(id(77) pos(3:-1:24565108 -> 3:-1:24566180)) length(1072)
         cluster(0) foldback be1(113: start 3:1:25401059) be2(119: start 3:1:25400602) length(457)
         cluster(3) foldback SV(id(88) pos(3:-1:26431918 -> 6:-1:26194040) : BND) with own breakend(88: start 3:-1:26431918) chain(0)
         cluster(4) foldback inversion SV(id(79) pos(3:1:26663922 -> 3:1:26664498)) length(576)
         cluster(0) SV(id(113) pos(3:1:25401059 -> 10:-1:60477224)) and cluster(1) SV(id(77) pos(3:-1:24565108 -> 3:-1:24566180)) have foldbacks on same arm
         cluster(0) SV(id(113) pos(3:1:25401059 -> 10:-1:60477224)) and cluster(3) SV(id(88) pos(3:-1:26431918 -> 6:-1:26194040)) have foldbacks on same arm
         cluster(0) SV(id(113) pos(3:1:25401059 -> 10:-1:60477224)) and cluster(4) SV(id(79) pos(3:1:26663922 -> 3:1:26664498)) have foldbacks on same arm
         cluster(2) boundary breakend(78: start 3:-1:25331584) ploidy TI match with cluster(0) breakend(119: start 3:1:25400602)
         cluster(2 svs=1) merges in other cluster(0 svs=7) and adopts new ID
         cluster(0) SV(id(77) pos(3:-1:24565108 -> 3:-1:24566180)) ploidy multiple(2, ploidy(4 vs min=2.0)
         cluster(0) SV(id(88) pos(3:-1:26431918 -> 6:-1:26194040)) ploidy multiple(2, ploidy(4 vs min=2.0)
         cluster(0) SV(id(79) pos(3:1:26663922 -> 3:1:26664498)) ploidy multiple(2, ploidy(4 vs min=2.0)
         cluster(0) starting chaining with assemblyLinks(4) svCount(8)
         index(0) method(ASSEMBLY) adding linked pair(119 12:72666892:end & 120 12:72667075:end ploidy=2.1) to new chain(0) ploidy(2.1 unc=0.5)
         pair(113 10:60477224:end & 120 10:60477422:start) links second breakend to chain(0) as only exhausted connection
         index(1) method(ASSEMBLY) adding linked pair(113 10:60477224:end & 120 10:60477422:start ploidy=2.0) to existing chain(0) ploidy(2.1 unc=0.4)
         assembly multi-sgl-conn pair(88 6:26194040:end & 89 6:26194117:start) ploidy(1.5): first(ploidy=3.5 links=2) second(ploidy=1.5 links=1)
         index(2) method(ASSEMBLY) adding linked pair(88 6:26194040:end & 89 6:26194117:start ploidy=1.5) to new chain(1) ploidy(1.5 unc=0.6)
         assembly multi-sgl-conn pair(88 6:26194040:end & 89 6:26194406:end) ploidy(1.8): first(ploidy=2.0 links=2) second(ploidy=1.5 links=1)
         pair(88 6:26194040:end & 89 6:26194406:end) links second breakend to chain(1) as only exhausted connection
         index(3) method(ASSEMBLY) adding linked pair(88 6:26194040:end & 89 6:26194406:end ploidy=1.8) to existing chain(1) ploidy(1.7 unc=0.5)
         created 2 partial chains from 4 assembly links
         cluster(0) chromosomes(4) AP totalSegments(4 valid=4)
         chain(0) adding chained foldback breakends(119: start 3:1:25400602 & 113: start 3:1:25401059 ploidy(2.1))
         chain(1) adding chained foldback breakends(88: start 3:-1:26431918 & 88: start 3:-1:26431918 ploidy(1.7))
         type-D: foldback breakend(79: start 3:1:26663922) ploidy(4.2) split by other foldback pair(88: start 3:-1:26431918 & 88: start 3:-1:26431918 ploidy(1.7))
         type-B: foldback(79: end 3:1:26664498) ploidy(4.2) matched with breakend(77: end 3:-1:24566180) ploidy(4.0)
         index(4) method(FOLDBACK) adding linked pair(79 3:26664498:end & 77 3:24566180:end ploidy=4.1) to new chain(2) ploidy(4.1 unc=0.4)
         end breakend exhausted: id(77) ploidy(3.5-4.0-4.5) counts(s=0.0 e=4.1)
         end breakend exhausted: id(79) ploidy(3.4-4.2-5.0) counts(s=0.0 e=4.1)
         foldback pair(79: start 3:1:26663922 & 79: end 3:1:26664498 ploidy(4.2)) removed from consideration
         foldback pair(77: start 3:-1:24565108 & 77: end 3:-1:24566180 ploidy(4.0)) removed from consideration
         type-B: foldback(113: start 3:1:25401059) ploidy(2.1) matched with breakend(78: start 3:-1:25331584) ploidy(2.0)
         type-A: foldback breakends(119: start 3:1:25400602 & 113: start 3:1:25401059) ploidy(2.1) exact split of breakend(77: start 3:-1:24565108) ploidy(4.1)
         type-A: foldback breakends(88: start 3:-1:26431918 & 88: start 3:-1:26431918) ploidy(1.7) exact split of breakend(79: start 3:1:26663922) ploidy(4.1)
         duplicating chain(2 links=1 sv=2) for multi-connect FOLDBACK_SPLIT
         index(5) method(FOLDBACK_SPLIT) adding linked pair(88 3:26431918:start & 79 3:26663922:start ploidy=1.9) to existing chain(2) ploidy(2.0 unc=0.3)
         index(5) method(FOLDBACK_SPLIT) adding linked pair(88 3:26431918:start & 79 3:26663922:start ploidy=1.9) to existing chain(2) ploidy(2.0 unc=0.3)
         start breakend exhausted: id(88) ploidy(2.6-3.5-4.4) counts(s=3.5 e=3.3)
         start breakend exhausted: id(79) ploidy(3.4-4.2-5.0) counts(s=4.2 e=4.1)
         SV(id(88) ploidy(2.6-3.5-4.4) counts(s=3.5 e=3.3)) both breakends exhausted
         SV(id(79) ploidy(3.4-4.2-5.0) counts(s=4.2 e=4.1)) both breakends exhausted
         foldback pair(88: start 3:-1:26431918 & 88: start 3:-1:26431918 ploidy(1.7)) removed from consideration
         chain(2) adding chained foldback breakends(77: start 3:-1:24565108 & 77: start 3:-1:24565108 ploidy(2.0))
         type-B: foldback(113: start 3:1:25401059) ploidy(2.1) matched with breakend(78: start 3:-1:25331584) ploidy(2.0)
         pair(113 3:25401059:start & 78 3:25331584 SGL-on-known) links first breakend to chain(0) as only exhausted connection
         index(6) method(FOLDBACK) adding linked pair(113 3:25401059:start & 78 3:25331584 SGL-on-known ploidy=2.0) to existing chain(0) ploidy(2.0 unc=0.4)
         start breakend exhausted: id(78) ploidy(1.0-2.0-2.9) counts(s=2.0 e=0.0)
         start breakend exhausted: id(113) ploidy(1.2-1.9-2.6) counts(s=1.9 e=2.1)
         SV(id(113) ploidy(1.2-1.9-2.6) counts(s=1.9 e=2.1)) both breakends exhausted
         foldback pair(119: start 3:1:25400602 & 113: start 3:1:25401059 ploidy(2.1)) removed from consideration
         type-B: foldback(77: start 3:-1:24565108) ploidy(2.0) matched with breakend(119: start 3:1:25400602) ploidy(2.0)
         pair(77 3:24565108:start & 119 3:25400602:start) links first breakend to chain(2) as only exhausted connection
         index(7) method(FOLDBACK) adding linked pair(77 3:24565108:start & 119 3:25400602:start ploidy=2.0) to existing chain(2) ploidy(2.0 unc=0.3)
         merging chain(0 links=3) start to chain(2 links=7) start
         start breakend exhausted: id(77) ploidy(3.5-4.0-4.5) counts(s=4.0 e=4.1)
         start breakend exhausted: id(119) ploidy(1.3-2.2-3.1) counts(s=2.2 e=2.1)
         SV(id(77) ploidy(3.5-4.0-4.5) counts(s=4.0 e=4.1)) both breakends exhausted
         SV(id(119) ploidy(1.3-2.2-3.1) counts(s=2.2 e=2.1)) both breakends exhausted
         foldback pair(77: start 3:-1:24565108 & 77: start 3:-1:24565108 ploidy(2.0)) removed from consideration
         cluster(0) chaining finished: chains(1 unique=1 links=8) SVs(8) unlinked SVs(0 ploidy=2.0) breakends(1 ploidy=2.0)

         cluster(0) added chain(0) ploidy(2.0) with 10 linked pairs:
         chain(0): 3_P_C - s_77_e - e_79_s - s_88_e - s_89_e - e_88_s - s_79_e - e_77_s - s_119_e - e_120_s - e_113_s - s_78_e - sgl_unclear
         chain(0) 0: pair(77 3:24566180:end & 79 3:26664498:end) FOLDBACK length(2098318) index(4)
         chain(0) 1: pair(79 3:26663922:start & 88 3:26431918:start) FOLDBACK_SPLIT length(232004) index(5)
         chain(0) 2: pair(88 6:26194040:end & 89 6:26194117:start) ASSEMBLY length(77) index(2)
         chain(0) 3: pair(89 6:26194406:end & 88 6:26194040:end) ASSEMBLY length(366) index(3)
         chain(0) 4: pair(88 3:26431918:start & 79 3:26663922:start) FOLDBACK_SPLIT length(232004) index(5)
         chain(0) 5: pair(79 3:26664498:end & 77 3:24566180:end) FOLDBACK length(2098318) index(4)
         chain(0) 6: pair(77 3:24565108:start & 119 3:25400602:start) FOLDBACK length(835494) index(7)
         chain(0) 7: pair(119 12:72666892:end & 120 12:72667075:end) ASSEMBLY length(183) index(0)
         chain(0) 8: pair(120 10:60477422:start & 113 10:60477224:end) ASSEMBLY length(198) index(1)
         chain(0) 9: pair(113 3:25401059:start & 78 3:25331584 SGL-on-known) FOLDBACK length(69475) index(6)
         cluster(0) complex SVs(8) desc(BND=4_INV=3_SGL=1 res=COMPLEX) arms(4) consis(1) chains(1 perc=1.00) replic(true) foldbacks=5 inv=3
         cluster(0) incomplete chains(1 incons=1) chainEnds(arms=0 repeats=0) unlinkedSVs(0 armCount(4 incons=0)) tiCount(short=0 long=0)        */
    }

    @Test
    public void testActualSimpleChaining1()
    {
        // based on a sample where the shortest TI of 76 bases is actually ignored so as to make 2 chains
        LinxTester tester = new LinxTester();

        final List<SvVarData> svList = SampleDataLoader.loadSampleTestData("CT_SAMPLE1");
        tester.AllVariants.addAll(svList);

        final SvVarData var1 = svList.get(0);
        final SvVarData var4 = svList.get(3);

        tester.CnDataLoader.getLohData().add(new LohEvent( "18", 23601785, 23577410,
                "BND", "BND", 1, var1.id(), var4.id()));

        tester.Analyser.getState().setSampleCnEventData(tester.CnDataLoader.getLohData(), tester.CnDataLoader.getHomLossData());
        tester.preClusteringInit();

        tester.Analyser.clusterAndAnalyse();

        // now check final chain-finding across all sub-clusters
        assertEquals(tester.Analyser.getClusters().size(), 1);
        final SvCluster cluster = tester.Analyser.getClusters().get(0);

        assertEquals(2, cluster.getChains().size());
    }

    @Test
    public void testActualFoldbackChaining1()
    {
        /*
         cluster(26) chaining finished: chains(1 unique=1 links=5) SVs(5) unlinked SVs(0 ploidy=2.0) breakends(2 ploidy=2.0)
         cluster(26) added chain(0) ploidy(1.0) with 5 linked pairs:
         chain(0): 18_P_T - s_128_e - e_129_s - e_127_s - s_129_e - s_125_e - s_130_e - 18_Q_C
         chain(0) 0: pair(128 18:20011830:end & 129 18:20306336:end) FOLDBACK length(294506) index(1)
         chain(0) 1: pair(129 18:20304431:start & 127 18:19821103:end) FOLDBACK_SPLIT length(483328) index(0)
         chain(0) 2: pair(127 18:19810967:start & 129 18:20304431:start) FOLDBACK_SPLIT length(493464) index(0)
         chain(0) 3: pair(129 18:20306336:end & 125 18:19019042:start) ONLY length(1287294) index(3)
         chain(0) 4: pair(125 18:19019748:end & 130 18:22311615:start) FOLDBACK length(3291867) index(2)

         */

        LinxTester tester = new LinxTester();

        final List<SvVarData> svList = SampleDataLoader.loadSampleTestData("FB_SAMPLE1");
        tester.AllVariants.addAll(svList);

        tester.preClusteringInit();
        tester.Analyser.clusterAndAnalyse();

        // now check final chain-finding across all sub-clusters
        assertEquals(tester.Analyser.getClusters().size(), 1);
        final SvCluster cluster = tester.Analyser.getClusters().get(0);

        assertEquals(4, cluster.getFoldbacks().size());
        assertTrue(cluster.getClusteringReasons().contains(FOLDBACKS.toString()));

        assertEquals(1, cluster.getChains().size());

        final SvChain chain = cluster.getChains().get(0);
        assertEquals(5, chain.getSvCount());
        assertEquals(5, chain.getLinkCount());
    }

    @Test
    public void testActualDoubleMinuteChaining()
    {
        // based on an EGFR AMP

        LinxTester tester = new LinxTester();

        // Configurator.setRootLevel(Level.DEBUG);

        /* SOLUTION: a closed loop:

             cluster(0) added chain(0) ploidy(2.0) with 8 linked pairs:
             chain(0): 7_P_T - s_83_e - e_84_s - s_80_e - e_81_s - e_79_s - e_78_s - s_82_e - e_84_s - s_83_e - 7_P_C
             chain(0) 0: pair(83 7:55636000:end & 84 7:55700000:end) ADJACENT length(64000) index(7)
             chain(0) 1: pair(84 7:54800000:start & 80 7:55200000:start) ONLY length(400000) index(6)
             chain(0) 2: pair(80 7:55223000:end & 81 7:55293000:end) PLOIDY_MATCH length(70000) index(5)
             chain(0) 3: pair(81 7:55204000:start & 79 7:55207000:end) PLOIDY_MATCH length(3000) index(3)
             chain(0) 4: pair(79 7:55145000:start & 78 7:55100000:end) PLOIDY_MATCH length(45000) index(4)
             chain(0) 5: pair(78 7:54877000:start & 82 7:55092000:start) ONLY length(215000) index(2)
             chain(0) 6: pair(82 7:55588000:end & 84 7:55700000:end) ONLY length(112000) index(1)
             chain(0) 7: pair(84 7:54800000:start & 83 7:54832000:start) ONLY length(32000) index(0)
         */

        final List<SvVarData> svList = SampleDataLoader.loadSampleTestData("DM_SAMPLE1");

        // final SvVarData dmDup = svList.get(6);

        tester.AllVariants.addAll(svList);

        tester.preClusteringInit();

        tester.Analyser.clusterAndAnalyse();

        // check clustering
        assertEquals(tester.Analyser.getClusters().size(), 1);

        final SvCluster cluster = tester.Analyser.getClusters().get(0);
        assertTrue(cluster.getAnnotations().contains(CLUSTER_ANNOT_DM));
        assertEquals(1, cluster.getDoubleMinuteSVs().size());

        // check chains
        assertEquals(2, cluster.getChains().size());
        final SvChain chain = cluster.getChains().get(0);
        assertEquals(8, chain.getLinkCount());
        assertEquals(7, chain.getSvCount());
        assertTrue(chain.isClosedLoop());
    }

    @Test
    public void testActualDoubleMinuteChaining2()
    {
        // based on an MDM2 AMP - ends up with 1 high-ploidy chain, another with some extra low-ploidy SVs, and 2 more disconnected

        LinxTester tester = new LinxTester();

        final List<SvVarData> svList = SampleDataLoader.loadSampleTestData("DM_SAMPLE2");

        tester.AllVariants.addAll(svList);

        tester.addLohEvent(svList.get(1).getBreakend(true), svList.get(13).getBreakend(true));

        tester.preClusteringInit();

        tester.Analyser.clusterAndAnalyse();

        // check clustering - second cluster is separate due to ploidy diffs
        assertEquals(2, tester.Analyser.getClusters().size());

        final SvCluster cluster = tester.Analyser.getClusters().get(0);
        assertTrue(cluster.getAnnotations().contains(CLUSTER_ANNOT_DM));
        assertEquals(6, cluster.getDoubleMinuteSVs().size());

        // check chains
        assertEquals(3, cluster.getChains().size());

        assertTrue(cluster.getUnlinkedSVs().isEmpty());

        final List<SvChain> dmChains = cluster.getDoubleMinuteChains();
        assertEquals(1, dmChains.size());
        final SvChain dmChain = dmChains.get(0);

        assertTrue(cluster.getChains().stream().anyMatch(x -> identicalChain(x, dmChain, false, true)));
   }
}
