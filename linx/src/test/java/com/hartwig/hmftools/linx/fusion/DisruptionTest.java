package com.hartwig.hmftools.linx.fusion;

import static com.hartwig.hmftools.common.test.GeneTestUtils.addGeneData;
import static com.hartwig.hmftools.common.test.GeneTestUtils.addTransExonData;
import static com.hartwig.hmftools.common.test.GeneTestUtils.createEnsemblGeneData;
import static com.hartwig.hmftools.common.test.GeneTestUtils.createGeneDataCache;
import static com.hartwig.hmftools.common.test.GeneTestUtils.createTransExons;
import static com.hartwig.hmftools.linx.utils.GeneTestUtils.CHR_1;
import static com.hartwig.hmftools.linx.utils.GeneTestUtils.CHR_2;
import static com.hartwig.hmftools.linx.utils.GeneTestUtils.GENE_ID_1;
import static com.hartwig.hmftools.linx.utils.GeneTestUtils.GENE_ID_2;
import static com.hartwig.hmftools.linx.utils.SampleDataLoader.setSvGeneData;
import static com.hartwig.hmftools.linx.utils.SvTestUtils.createBnd;
import static com.hartwig.hmftools.linx.utils.SvTestUtils.createDel;
import static com.hartwig.hmftools.linx.utils.SvTestUtils.createDup;
import static com.hartwig.hmftools.linx.utils.SvTestUtils.createIns;
import static com.hartwig.hmftools.linx.utils.SvTestUtils.createInv;
import static com.hartwig.hmftools.linx.utils.SvTestUtils.createSgl;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertFalse;
import static org.junit.Assert.assertTrue;

import java.util.List;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.ensemblcache.EnsemblDataCache;
import com.hartwig.hmftools.common.gene.GeneData;
import com.hartwig.hmftools.common.gene.TranscriptData;
import com.hartwig.hmftools.linx.types.SvVarData;
import com.hartwig.hmftools.linx.utils.LinxTester;

import org.junit.Test;

public class DisruptionTest
{
    @Test
    public void testSimpleDisruptions()
    {
        LinxTester tester = new LinxTester();

        EnsemblDataCache geneTransCache = createGeneDataCache();
        tester.initialiseFusions(geneTransCache);

        String chromosome = CHR_1;
        String geneId = GENE_ID_1;
        addTestGeneData(geneTransCache, chromosome, geneId);

        // non-disruptive simple SVs
        SvVarData var1 = createDel(tester.nextVarId(), chromosome, 11000,12000);
        SvVarData var2 = createDup(tester.nextVarId(), chromosome, 21000,22000);
        SvVarData var3 = createIns(tester.nextVarId(), chromosome, 31000,32000);

        // disruptive other SVs
        SvVarData var4 = createInv(tester.nextVarId(), chromosome, 41000, 42000, -1);
        SvVarData var5 = createBnd(tester.nextVarId(), chromosome, 51000, -1, "2", 1000, -1);
        SvVarData var6 = createSgl(tester.nextVarId(), chromosome, 61000, -1);

        // DELs and DUPs crossing exon boundaries
        SvVarData var7 = createDel(tester.nextVarId(), chromosome, 68000,72000);
        SvVarData var8 = createDup(tester.nextVarId(), chromosome, 78000,96000);

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

        // assertEquals(8, tester.Analyser.getClusters().size());

        final DisruptionFinder disruptionFinder = tester.FusionAnalyser.getDisruptionFinder();
        disruptionFinder.addDisruptionGene(geneId);

        setSvGeneData(tester.AllVariants, geneTransCache, true);
        tester.FusionAnalyser.annotateTranscripts(tester.AllVariants, true);

        tester.AllVariants.forEach(x -> assertEquals(1, x.getGenesList(true).size()));

        assertFalse(var1.getGenesList(true).get(0).transcripts().get(0).isDisruptive());
        assertFalse(var2.getGenesList(true).get(0).transcripts().get(0).isDisruptive());
        assertFalse(var3.getGenesList(true).get(0).transcripts().get(0).isDisruptive());

        assertTrue(var4.getGenesList(true).get(0).transcripts().get(0).isDisruptive());
        assertTrue(var5.getGenesList(true).get(0).transcripts().get(0).isDisruptive());
        assertTrue(var6.getGenesList(true).get(0).transcripts().get(0).isDisruptive());
        assertTrue(var7.getGenesList(true).get(0).transcripts().get(0).isDisruptive());
        assertTrue(var8.getGenesList(true).get(0).transcripts().get(0).isDisruptive());

        disruptionFinder.findReportableDisruptions(tester.AllVariants, tester.Analyser.getClusters());
        assertEquals(8, disruptionFinder.getDisruptions().size());

        tester.clearClustersAndSVs();

        // a DUP repeating the first exon (which has no splice acceptor) is not disruptive - unless it is chained
        SvVarData var9 = createDup(tester.nextVarId(), chromosome, 10050,10500);

        tester.AllVariants.add(var9);

        tester.preClusteringInit();
        tester.Analyser.clusterAndAnalyse();
        setSvGeneData(tester.AllVariants, geneTransCache, false);
        tester.FusionAnalyser.annotateTranscripts(tester.AllVariants, true);

        assertFalse(var9.getGenesList(true).get(0).transcripts().get(0).isDisruptive());
        assertFalse(var9.getGenesList(false).get(0).transcripts().get(0).isDisruptive());

        disruptionFinder.findReportableDisruptions(tester.AllVariants, tester.Analyser.getClusters());
        assertTrue(disruptionFinder.getDisruptions().isEmpty());

        tester.clearClustersAndSVs();

        SvVarData var10 = createDup(tester.nextVarId(), chromosome, 9500,10500);
        tester.AllVariants.add(var10);

        tester.preClusteringInit();
        tester.Analyser.clusterAndAnalyse();
        setSvGeneData(tester.AllVariants, geneTransCache, false);
        tester.FusionAnalyser.annotateTranscripts(tester.AllVariants, true);

        assertTrue(var10.getGenesList(true).isEmpty());
        assertFalse(var10.getGenesList(false).get(0).transcripts().get(0).isDisruptive());

        disruptionFinder.findReportableDisruptions(tester.AllVariants, tester.Analyser.getClusters());
        assertTrue(disruptionFinder.getDisruptions().isEmpty());

        tester.clearClustersAndSVs();

        // DUPs which have only 1 end in a transcript or gene
        SvVarData var11 = createDup(tester.nextVarId(), chromosome, 200,10500);
        tester.AllVariants.add(var11);

        tester.preClusteringInit();
        tester.Analyser.clusterAndAnalyse();
        setSvGeneData(tester.AllVariants, geneTransCache, false);
        tester.FusionAnalyser.annotateTranscripts(tester.AllVariants, true);

        assertFalse(var11.getGenesList(false).get(0).transcripts().get(0).isDisruptive());

        disruptionFinder.findReportableDisruptions(tester.AllVariants, tester.Analyser.getClusters());
        assertTrue(disruptionFinder.getDisruptions().isEmpty());

        tester.clearClustersAndSVs();

        SvVarData var12 = createDup(tester.nextVarId(), chromosome, 10500,15000);
        tester.AllVariants.add(var12);

        tester.preClusteringInit();
        tester.Analyser.clusterAndAnalyse();
        setSvGeneData(tester.AllVariants, geneTransCache, false);
        tester.FusionAnalyser.annotateTranscripts(tester.AllVariants, true);

        assertFalse(var12.getGenesList(true).get(0).transcripts().get(0).isDisruptive());

        disruptionFinder.findReportableDisruptions(tester.AllVariants, tester.Analyser.getClusters());
        assertTrue(disruptionFinder.getDisruptions().isEmpty());
    }

    private void addTestGeneData(EnsemblDataCache geneTransCache, final String chromosome, final String geneId)
    {
        byte strand = 1;

        List<GeneData> geneList = Lists.newArrayList();
        geneList.add(createEnsemblGeneData(geneId, geneId, chromosome, strand, 10000, 120000));
        addGeneData(geneTransCache, chromosome, geneList);

        List<TranscriptData> transDataList = Lists.newArrayList();

        int transId = 1;

        int[] exonStarts = new int[] { 10000, 20000, 30000, 40000, 50000, 60000, 70000, 80000, 90000, 100000, 110000 };

        int codingStart = 20002;
        int codingEnd = 110098;
        TranscriptData transData = createTransExons(geneId, transId++, strand, exonStarts,  100, codingStart, codingEnd, true, "");
        transDataList.add(transData);

        addTransExonData(geneTransCache, geneId, transDataList);
    }

    @Test
    public void testChainedDisruptions()
    {
        LinxTester tester = new LinxTester();

        EnsemblDataCache geneTransCache = createGeneDataCache();
        tester.initialiseFusions(geneTransCache);

        String chromosome = CHR_1;
        String geneId = GENE_ID_1;

        String chromosome2 = CHR_2;
        String geneId2 = GENE_ID_2;

        addTestGeneData(geneTransCache, chromosome, geneId);
        addTestGeneData(geneTransCache, chromosome2, geneId2);

        // scenarios (TRUE - disruptive, FALSE - not disruptive)
        // - FALSE: chain which goes out and back without traversing another splice acceptor, ends on correct orientation
        // - TRUE: chain which goes out and back but traverses another splice acceptor
        // - TRUE: chain which goes out and back without traversing another splice acceptor, ends on incorrect orientation

        // test 1: a pair of BNDs forming a remote TI contained within an intron
        SvVarData var1 = createBnd(tester.nextVarId(), chromosome, 11000, 1, "3", 100, -1);
        SvVarData var2 = createBnd(tester.nextVarId(), chromosome, 12000, -1, "3", 200, 1);

        tester.AllVariants.add(var1);
        tester.AllVariants.add(var2);

        prepareGeneAnnotations(tester, geneTransCache);

        assertEquals(1, tester.Analyser.getClusters().size());

        tester.AllVariants.forEach(x -> assertEquals(1, x.getGenesList(true).size()));

        assertTrue(!var1.getGenesList(true).get(0).transcripts().get(0).isDisruptive());
        assertTrue(!var2.getGenesList(true).get(0).transcripts().get(0).isDisruptive());


        // test 2: same again but with 2 INVs (reciprocal INV) within an intron
        tester.clearClustersAndSVs();

        var1 = createInv(tester.nextVarId(), chromosome, 11000, 14000, 1);
        var2 = createInv(tester.nextVarId(), chromosome, 13900, 18000, -1);

        tester.AllVariants.add(var1);
        tester.AllVariants.add(var2);

        prepareGeneAnnotations(tester, geneTransCache);

        tester.AllVariants.forEach(x -> assertEquals(1, x.getGenesList(true).size()));

        assertTrue(!var1.getGenesList(true).get(0).transcripts().get(0).isDisruptive());
        assertTrue(!var1.getGenesList(false).get(0).transcripts().get(0).isDisruptive());
        assertTrue(!var2.getGenesList(true).get(0).transcripts().get(0).isDisruptive());
        assertTrue(!var2.getGenesList(false).get(0).transcripts().get(0).isDisruptive());


        // test 3: same again but with a longer chain
        tester.clearClustersAndSVs();

        var1 = createDel(tester.nextVarId(), chromosome, 11000, 12000);
        var2 = createBnd(tester.nextVarId(), chromosome, 13000, 1, "3", 100, -1);
        SvVarData var3 = createBnd(tester.nextVarId(), chromosome, 14000, -1, "3", 200, 1);
        SvVarData var4 = createDel(tester.nextVarId(), chromosome, 15000, 16000);

        tester.AllVariants.add(var1);
        tester.AllVariants.add(var2);
        tester.AllVariants.add(var3);
        tester.AllVariants.add(var4);

        prepareGeneAnnotations(tester, geneTransCache);

        tester.AllVariants.forEach(x -> assertEquals(1, x.getGenesList(true).size()));

        assertTrue(!var1.getGenesList(true).get(0).transcripts().get(0).isDisruptive());
        assertTrue(!var2.getGenesList(true).get(0).transcripts().get(0).isDisruptive());
        assertTrue(!var3.getGenesList(true).get(0).transcripts().get(0).isDisruptive());
        assertTrue(!var4.getGenesList(true).get(0).transcripts().get(0).isDisruptive());


        // test 4: now with invalid final orientation - even though it makes an intronic TI
        tester.clearClustersAndSVs();

        var1 = createDel(tester.nextVarId(), chromosome, 11000, 12000);
        var2 = createBnd(tester.nextVarId(), chromosome, 13000, 1, "3", 100, -1);
        var3 = createBnd(tester.nextVarId(), chromosome, 14000, 1, "3", 200, 1);
        var4 = createDup(tester.nextVarId(), chromosome, 13500, 16000);

        tester.AllVariants.add(var1);
        tester.AllVariants.add(var2);
        tester.AllVariants.add(var3);
        tester.AllVariants.add(var4);

        prepareGeneAnnotations(tester, geneTransCache);

        tester.AllVariants.forEach(x -> assertEquals(1, x.getGenesList(true).size()));

        assertTrue(!var1.getGenesList(true).get(0).transcripts().get(0).isDisruptive());
        assertTrue(!var1.getGenesList(false).get(0).transcripts().get(0).isDisruptive());
        assertTrue(var2.getGenesList(true).get(0).transcripts().get(0).isDisruptive());
        assertTrue(var3.getGenesList(true).get(0).transcripts().get(0).isDisruptive());
        assertTrue(!var4.getGenesList(true).get(0).transcripts().get(0).isDisruptive());


        // now with an invalid traversal
        tester.clearClustersAndSVs();

        var1 = createDel(tester.nextVarId(), chromosome, 11000, 12000);
        var2 = createBnd(tester.nextVarId(), chromosome, 13000, 1, chromosome2, 18000, -1);
        var3 = createBnd(tester.nextVarId(), chromosome2, 22000, 1, "3", 1000, -1);
        var4 = createBnd(tester.nextVarId(), chromosome, 14000, -1, "3", 3000, 1);
        SvVarData var5 = createDel(tester.nextVarId(), chromosome, 15000, 16000);

        tester.AllVariants.add(var1);
        tester.AllVariants.add(var2);
        tester.AllVariants.add(var3);
        tester.AllVariants.add(var4);
        tester.AllVariants.add(var5);

        prepareGeneAnnotations(tester, geneTransCache);

        assertEquals(1, tester.Analyser.getClusters().size());

        tester.AllVariants.forEach(x -> assertEquals(1, x.getGenesList(true).size()));

        assertTrue(!var1.getGenesList(true).get(0).transcripts().get(0).isDisruptive());

        assertTrue(var2.getGenesList(true).get(0).transcripts().get(0).isDisruptive());
        assertTrue(var2.getGenesList(false).get(0).transcripts().get(0).isDisruptive());

        assertTrue(var3.getGenesList(true).get(0).transcripts().get(0).isDisruptive());

        assertTrue(var4.getGenesList(true).get(0).transcripts().get(0).isDisruptive());

        assertTrue(!var5.getGenesList(true).get(0).transcripts().get(0).isDisruptive());
    }

    @Test
    public void testChainedDisruptions2()
    {
        LinxTester tester = new LinxTester();

        EnsemblDataCache geneTransCache = createGeneDataCache();
        tester.initialiseFusions(geneTransCache);

        String chromosome = CHR_1;
        String geneId = GENE_ID_1;

        addTestGeneData(geneTransCache, chromosome, geneId);

        String chromosome2 = CHR_2;
        String geneId2 = GENE_ID_2;

        addTestGeneData(geneTransCache, chromosome2, geneId2);

        // test 1: pair of INVs forming a remote TI contained within an intron, with the other ends non-genic
        SvVarData var1 = createInv(tester.nextVarId(), chromosome, 1000, 15000, 1);
        SvVarData var2 = createInv(tester.nextVarId(), chromosome, 2000, 14000, -1);

        tester.AllVariants.add(var1);
        tester.AllVariants.add(var2);

        prepareGeneAnnotations(tester, geneTransCache);

        assertEquals(1, tester.Analyser.getClusters().size());

        tester.AllVariants.forEach(x -> assertEquals(1, x.getGenesList(false).size()));

        assertTrue(!var1.getGenesList(false).get(0).transcripts().get(0).isDisruptive());
        assertTrue(!var2.getGenesList(false).get(0).transcripts().get(0).isDisruptive());

        tester.clearClustersAndSVs();


        // test 2: now with one of the variants ending in another gene - makes the first INV also disruptive
        var1 = createInv(tester.nextVarId(), chromosome, 1000, 15000, 1);
        var2 = createBnd(tester.nextVarId(), chromosome, 14000, -1, chromosome2, 25000, -1);

        tester.AllVariants.add(var1);
        tester.AllVariants.add(var2);

        prepareGeneAnnotations(tester, geneTransCache);

        assertEquals(1, tester.Analyser.getClusters().size());

        tester.AllVariants.forEach(x -> assertEquals(1, x.getGenesList(false).size()));

        assertTrue(var1.getGenesList(true).isEmpty());
        assertTrue(var1.getGenesList(false).get(0).transcripts().get(0).isDisruptive());
        assertTrue(var2.getGenesList(true).get(0).transcripts().get(0).isDisruptive());
        assertTrue(var2.getGenesList(false).get(0).transcripts().get(0).isDisruptive());


        // test 3: small shattering within an intron and a remote TI - non disruptive
        tester.clearClustersAndSVs();

        String remoteChr = "3";
        var1 = createInv(tester.nextVarId(), chromosome, 12000, 14000, 1);
        var2 = createBnd(tester.nextVarId(), chromosome, 13000, -1, remoteChr, 100, -1);

        SvVarData var3 = createBnd(tester.nextVarId(), chromosome, 17000, 1, remoteChr, 200, 1);
        SvVarData var4 = createInv(tester.nextVarId(), chromosome, 16000, 18000, -1);

        tester.AllVariants.add(var1);
        tester.AllVariants.add(var2);
        tester.AllVariants.add(var3);
        tester.AllVariants.add(var4);

        prepareGeneAnnotations(tester, geneTransCache);

        assertEquals(1, tester.Analyser.getClusters().size());

        assertTrue(!var1.getGenesList(true).get(0).transcripts().get(0).isDisruptive());
        assertTrue(!var1.getGenesList(false).get(0).transcripts().get(0).isDisruptive());
        assertTrue(!var2.getGenesList(true).get(0).transcripts().get(0).isDisruptive());
        assertTrue(!var3.getGenesList(true).get(0).transcripts().get(0).isDisruptive());
        assertTrue(!var4.getGenesList(true).get(0).transcripts().get(0).isDisruptive());


        // test 4: multiple TIs within the same intro within a chain which has its ends outside genic regions
        // sorounded by a non-disruptive intronic DEL
        tester.clearClustersAndSVs();

        var1 = createDel(tester.nextVarId(), chromosome, 11800, 14500);
        var2 = createBnd(tester.nextVarId(), chromosome, 12000, -1, remoteChr, 100, 1);
        var3 = createInv(tester.nextVarId(), chromosome, 12100, 14100, 1);

        var2.setAssemblyData(true, "asmb_A2_A3");
        var3.setAssemblyData(true, "asmb_A2_A3");

        var4 = createBnd(tester.nextVarId(), chromosome, 14000, -1, remoteChr, 200, -1);
        var3.setAssemblyData(false, "asmb_A3_A4");
        var4.setAssemblyData(true, "asmb_A3_A4");

        tester.AllVariants.add(var1);
        tester.AllVariants.add(var2);
        tester.AllVariants.add(var3);
        tester.AllVariants.add(var4);

        prepareGeneAnnotations(tester, geneTransCache);

        assertEquals(1, tester.Analyser.getClusters().size());

        assertTrue(!var1.getGenesList(true).get(0).transcripts().get(0).isDisruptive());
        assertTrue(!var1.getGenesList(false).get(0).transcripts().get(0).isDisruptive());

        assertTrue(!var3.getGenesList(true).get(0).transcripts().get(0).isDisruptive());
        assertTrue(!var3.getGenesList(false).get(0).transcripts().get(0).isDisruptive());
        assertTrue(!var4.getGenesList(true).get(0).transcripts().get(0).isDisruptive());
        assertTrue(!var3.getGenesList(true).get(0).transcripts().get(0).isDisruptive());


        // test 5: single remote TI without an intron, but forming DBs with a DEL in the same cluster
        // meaning the TI is not isolated and remain disruptive

        tester.clearClustersAndSVs();

        var1 = createDel(tester.nextVarId(), chromosome, 12000, 25500);
        var2 = createBnd(tester.nextVarId(), chromosome, 12100, -1, remoteChr, 100, 1);
        var3 = createBnd(tester.nextVarId(), chromosome, 12500, 1, remoteChr, 200, -1);

        var2.setAssemblyData(true, "asmb_A2_A3");
        var3.setAssemblyData(true, "asmb_A2_A3");

        tester.AllVariants.add(var1);
        tester.AllVariants.add(var2);
        tester.AllVariants.add(var3);

        prepareGeneAnnotations(tester, geneTransCache);

        assertEquals(1, tester.Analyser.getClusters().size());

        assertTrue(var1.getGenesList(true).get(0).transcripts().get(0).isDisruptive());
        assertTrue(var1.getGenesList(false).get(0).transcripts().get(0).isDisruptive());
        assertTrue(var2.getGenesList(true).get(0).transcripts().get(0).isDisruptive());
        assertTrue(var3.getGenesList(true).get(0).transcripts().get(0).isDisruptive());


        // test 5: more complicated example, with 2 TIs spanning an intron, but forming DBs with SVs in another chain in the same cluster
        // menaning the TIs are not isolated and remain disruptive

        tester.clearClustersAndSVs();

        var1 = createDel(tester.nextVarId(), chromosome, 12000, 25500);
        var2 = createBnd(tester.nextVarId(), chromosome, 12100, -1, remoteChr, 100, 1);
        var3 = createInv(tester.nextVarId(), chromosome, 12500, 25400, 1);
        var4 = createBnd(tester.nextVarId(), chromosome, 25000, -1, remoteChr, 200, -1);

        tester.AllVariants.add(var1);
        tester.AllVariants.add(var2);
        tester.AllVariants.add(var3);
        tester.AllVariants.add(var4);

        prepareGeneAnnotations(tester, geneTransCache);

        assertEquals(1, tester.Analyser.getClusters().size());

        assertTrue(var1.getGenesList(true).get(0).transcripts().get(0).isDisruptive());
        assertTrue(var1.getGenesList(false).get(0).transcripts().get(0).isDisruptive());
        assertTrue(var2.getGenesList(true).get(0).transcripts().get(0).isDisruptive());
        assertTrue(var3.getGenesList(true).get(0).transcripts().get(0).isDisruptive());
        assertTrue(var3.getGenesList(false).get(0).transcripts().get(0).isDisruptive());
        assertTrue(var4.getGenesList(true).get(0).transcripts().get(0).isDisruptive());

    }

    private void prepareGeneAnnotations(final LinxTester tester, final EnsemblDataCache geneTransCache)
    {
        tester.preClusteringInit();
        tester.Analyser.clusterAndAnalyse();
        setSvGeneData(tester.AllVariants, geneTransCache, false);
        tester.FusionAnalyser.annotateTranscripts(tester.AllVariants, true);
    }
}
