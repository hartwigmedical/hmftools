package com.hartwig.hmftools.linx.fusion;

import static com.hartwig.hmftools.common.test.GeneTestUtils.addGeneData;
import static com.hartwig.hmftools.common.test.GeneTestUtils.addTransExonData;
import static com.hartwig.hmftools.common.test.GeneTestUtils.createEnsemblGeneData;
import static com.hartwig.hmftools.common.test.GeneTestUtils.createGeneDataCache;
import static com.hartwig.hmftools.common.test.GeneTestUtils.createTransExons;
import static com.hartwig.hmftools.linx.analysis.VariantPrep.setSvGeneData;
import static com.hartwig.hmftools.linx.utils.GeneTestUtils.CHR_1;
import static com.hartwig.hmftools.linx.utils.GeneTestUtils.CHR_2;
import static com.hartwig.hmftools.linx.utils.GeneTestUtils.GENE_ID_1;
import static com.hartwig.hmftools.linx.utils.GeneTestUtils.GENE_ID_2;
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
import com.hartwig.hmftools.linx.SampleAnalyser;
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

        // a DUP repeating the first exon (which has no splice acceptor)
        SvVarData var9 = createDup(tester.nextVarId(), chromosome, 10050,10500);
        SvVarData var10 = createDup(tester.nextVarId(), chromosome, 9500,10500);

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

        tester.preClusteringInit();
        tester.Analyser.clusterAndAnalyse();

        assertEquals(8, tester.Analyser.getClusters().size());

        final DisruptionFinder disruptionFinder = tester.FusionAnalyser.getDisruptionFinder();
        disruptionFinder.addDisruptionGene(geneTransCache.getGeneDataById(geneId));

        setSvGeneData(tester.AllVariants, geneTransCache, true, false);
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
        assertFalse(var9.getGenesList(true).get(0).transcripts().get(0).isDisruptive());
        assertFalse(var10.getGenesList(true).get(0).transcripts().get(0).isDisruptive());

        disruptionFinder.findReportableDisruptions(tester.AllVariants);
        assertEquals(8, disruptionFinder.getDisruptions().size());
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

        final DisruptionFinder disruptionFinder = tester.FusionAnalyser.getDisruptionFinder();

        // scenarios (TRUE - disruptive, FALSE - not disruptive)
        // - FALSE: chain which goes out and back without traversing another splice acceptor, ends on correct orientation
        // - TRUE: chain which goes out and back but traverses another splice acceptor
        // - TRUE: chain which goes out and back without traversing another splice acceptor, ends on incorrect orientation

        SvVarData var1 = createBnd(tester.nextVarId(), chromosome, 11000, 1, "3", 100, -1);
        SvVarData var2 = createBnd(tester.nextVarId(), chromosome, 12000, -1, "3", 200, 1);

        tester.AllVariants.add(var1);
        tester.AllVariants.add(var2);

        tester.preClusteringInit();
        tester.Analyser.clusterAndAnalyse();

        assertEquals(1, tester.Analyser.getClusters().size());

        setSvGeneData(tester.AllVariants, geneTransCache, false, false);
        tester.FusionAnalyser.annotateTranscripts(tester.AllVariants, true);

        tester.AllVariants.forEach(x -> assertEquals(1, x.getGenesList(true).size()));

        assertTrue(!var1.getGenesList(true).get(0).transcripts().get(0).isDisruptive());
        assertTrue(!var2.getGenesList(true).get(0).transcripts().get(0).isDisruptive());


        // same again but with a longer chain
        tester.clearClustersAndSVs();

        var1 = createDel(tester.nextVarId(), chromosome, 11000, 12000);
        var2 = createBnd(tester.nextVarId(), chromosome, 13000, 1, "3", 100, -1);
        SvVarData var3 = createBnd(tester.nextVarId(), chromosome, 14000, -1, "3", 200, 1);
        SvVarData var4 = createDel(tester.nextVarId(), chromosome, 15000, 16000);

        tester.AllVariants.add(var1);
        tester.AllVariants.add(var2);
        tester.AllVariants.add(var3);
        tester.AllVariants.add(var4);

        tester.preClusteringInit();
        tester.Analyser.clusterAndAnalyse();

        assertEquals(1, tester.Analyser.getClusters().size());

        setSvGeneData(tester.AllVariants, geneTransCache, false, false);
        tester.FusionAnalyser.annotateTranscripts(tester.AllVariants, true);

        tester.AllVariants.forEach(x -> assertEquals(1, x.getGenesList(true).size()));

        assertTrue(!var1.getGenesList(true).get(0).transcripts().get(0).isDisruptive());
        assertTrue(!var2.getGenesList(true).get(0).transcripts().get(0).isDisruptive());
        assertTrue(!var3.getGenesList(true).get(0).transcripts().get(0).isDisruptive());
        assertTrue(!var4.getGenesList(true).get(0).transcripts().get(0).isDisruptive());

        // now with invalid final orientation - even though it makes an intronic TI
        tester.clearClustersAndSVs();

        var1 = createDel(tester.nextVarId(), chromosome, 11000, 12000);
        var2 = createBnd(tester.nextVarId(), chromosome, 13000, 1, "3", 100, -1);
        var3 = createBnd(tester.nextVarId(), chromosome, 14000, 1, "3", 200, 1);
        var4 = createDup(tester.nextVarId(), chromosome, 13500, 16000);

        tester.AllVariants.add(var1);
        tester.AllVariants.add(var2);
        tester.AllVariants.add(var3);
        tester.AllVariants.add(var4);

        tester.preClusteringInit();
        tester.Analyser.clusterAndAnalyse();

        assertEquals(1, tester.Analyser.getClusters().size());

        setSvGeneData(tester.AllVariants, geneTransCache, false, false);
        tester.FusionAnalyser.annotateTranscripts(tester.AllVariants, true);

        tester.AllVariants.forEach(x -> assertEquals(1, x.getGenesList(true).size()));

        assertTrue(!var1.getGenesList(true).get(0).transcripts().get(0).isDisruptive());
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

        tester.preClusteringInit();
        tester.Analyser.clusterAndAnalyse();

        assertEquals(1, tester.Analyser.getClusters().size());

        setSvGeneData(tester.AllVariants, geneTransCache, false, false);
        tester.FusionAnalyser.annotateTranscripts(tester.AllVariants, true);

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

        final DisruptionFinder disruptionFinder = tester.FusionAnalyser.getDisruptionFinder();

        // scenarios (TRUE - disruptive, FALSE - not disruptive)
        // - FALSE: chain which starts outside genic regions and only has short TI inside same intron
        // - TRUE: chain which starts inside genic regions and only has short TI inside same intron

        SvVarData var1 = createInv(tester.nextVarId(), chromosome, 1000, 15000, 1);
        SvVarData var2 = createInv(tester.nextVarId(), chromosome, 2000, 14000, -1);

        tester.AllVariants.add(var1);
        tester.AllVariants.add(var2);

        tester.preClusteringInit();
        tester.Analyser.clusterAndAnalyse();

        assertEquals(1, tester.Analyser.getClusters().size());

        setSvGeneData(tester.AllVariants, geneTransCache, false, false);
        tester.FusionAnalyser.annotateTranscripts(tester.AllVariants, true);

        tester.AllVariants.forEach(x -> assertEquals(1, x.getGenesList(false).size()));

        assertTrue(!var1.getGenesList(false).get(0).transcripts().get(0).isDisruptive());
        assertTrue(!var2.getGenesList(false).get(0).transcripts().get(0).isDisruptive());

        tester.clearClustersAndSVs();

        // now with one of the variants ending in another gene
        var1 = createInv(tester.nextVarId(), chromosome, 1000, 15000, 1);
        var2 = createBnd(tester.nextVarId(), chromosome, 14000, -1, chromosome2, 25000, -1);

        tester.AllVariants.add(var1);
        tester.AllVariants.add(var2);

        tester.preClusteringInit();
        tester.Analyser.clusterAndAnalyse();

        assertEquals(1, tester.Analyser.getClusters().size());

        setSvGeneData(tester.AllVariants, geneTransCache, false, false);
        tester.FusionAnalyser.annotateTranscripts(tester.AllVariants, true);

        tester.AllVariants.forEach(x -> assertEquals(1, x.getGenesList(false).size()));

        assertTrue(var1.getGenesList(false).get(0).transcripts().get(0).isDisruptive());
        assertTrue(var2.getGenesList(false).get(0).transcripts().get(0).isDisruptive());
    }

}
