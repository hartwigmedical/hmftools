package com.hartwig.hmftools.linx.germline;

import static com.hartwig.hmftools.common.fusion.FusionCommon.POS_STRAND;
import static com.hartwig.hmftools.common.test.GeneTestUtils.addGeneData;
import static com.hartwig.hmftools.common.test.GeneTestUtils.addTransExonData;
import static com.hartwig.hmftools.common.test.GeneTestUtils.createEnsemblGeneData;
import static com.hartwig.hmftools.common.test.GeneTestUtils.createGeneDataCache;
import static com.hartwig.hmftools.common.test.GeneTestUtils.createTransExons;
import static com.hartwig.hmftools.common.utils.sv.SvCommonUtils.NEG_ORIENT;
import static com.hartwig.hmftools.common.utils.sv.SvCommonUtils.POS_ORIENT;
import static com.hartwig.hmftools.linx.gene.BreakendGenePrep.setSvGeneData;
import static com.hartwig.hmftools.linx.utils.GeneTestUtils.CHR_1;
import static com.hartwig.hmftools.linx.utils.GeneTestUtils.CHR_2;
import static com.hartwig.hmftools.linx.utils.GeneTestUtils.GENE_ID_1;
import static com.hartwig.hmftools.linx.utils.GeneTestUtils.GENE_ID_2;
import static com.hartwig.hmftools.linx.utils.GeneTestUtils.GENE_NAME_1;
import static com.hartwig.hmftools.linx.utils.GeneTestUtils.createDriverGene;
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
import com.hartwig.hmftools.common.drivercatalog.DriverCatalog;
import com.hartwig.hmftools.common.ensemblcache.EnsemblDataCache;
import com.hartwig.hmftools.common.gene.GeneData;
import com.hartwig.hmftools.common.gene.TranscriptData;
import com.hartwig.hmftools.common.sv.linx.LinxGermlineSv;
import com.hartwig.hmftools.linx.fusion.DisruptionFinder;
import com.hartwig.hmftools.linx.fusion.SvDisruptionData;
import com.hartwig.hmftools.linx.types.SglMapping;
import com.hartwig.hmftools.linx.types.SvVarData;
import com.hartwig.hmftools.linx.utils.LinxTester;

import org.junit.Test;

public class GermlineDisruptionTest
{
    private final LinxTester mLinx;

    private final DisruptionFinder mDisruptionFinder;
    private final GermlineDisruptions mGermlineDisruptions;

    public GermlineDisruptionTest()
    {
        mLinx = new LinxTester(true);

        mLinx.Config.DriverGenes.add(createDriverGene(GENE_NAME_1));

        EnsemblDataCache geneTransCache = createGeneDataCache();

        String chromosome = CHR_1;
        String geneId = GENE_ID_1;
        addTestGeneData(geneTransCache, chromosome, geneId, GENE_NAME_1);

        mLinx.initialiseFusions(geneTransCache);

        mDisruptionFinder = mLinx.FusionAnalyser.getDisruptionFinder();
        mDisruptionFinder.addDisruptionGene(geneId);

        mGermlineDisruptions = mDisruptionFinder.germlineDisruptions();
    }

    @Test
    public void testWholeGeneDeletions()
    {
        // DEL deleting whole gene
        SvVarData var1 = createDel(mLinx.nextVarId(), CHR_1, 5000,150000);

        mLinx.AllVariants.add(var1);

        mLinx.preClusteringInit();
        mLinx.Analyser.clusterAndAnalyse();

        mDisruptionFinder.findReportableDisruptions(mLinx.AllVariants, mLinx.Analyser.getClusters());

        List<LinxGermlineSv> germlineSVs = getGermlineSVs();
        assertEquals(1, germlineSVs.size());
        assertTrue(germlineSVs.get(0).Reported);
    }

    @Test
    public void testSglMappedDel()
    {
        // SGL with mapping to make a DEL around a gene
        SvVarData var = createSgl(mLinx.nextVarId(), CHR_1, 5000,POS_ORIENT);
        var.getSglMappings().add(new SglMapping(CHR_1, 150000, NEG_ORIENT, "", 1));

        mLinx.AllVariants.add(var);

        mLinx.preClusteringInit();
        mLinx.Analyser.clusterAndAnalyse();

        mDisruptionFinder.findReportableDisruptions(mLinx.AllVariants, mLinx.Analyser.getClusters());

        List<LinxGermlineSv> germlineSVs = getGermlineSVs();
        assertEquals(1, germlineSVs.size());
        assertTrue(germlineSVs.get(0).Reported);

        mLinx.clearClustersAndSVs();

        // partial DEL - exonic to intronic
        var = createSgl(mLinx.nextVarId(), CHR_1, 5000,POS_ORIENT);
        var.getSglMappings().add(new SglMapping(CHR_1, 12100, NEG_ORIENT, "", 1));

        mLinx.AllVariants.add(var);

        mLinx.preClusteringInit();
        mLinx.Analyser.clusterAndAnalyse();

        mDisruptionFinder.findReportableDisruptions(mLinx.AllVariants, mLinx.Analyser.getClusters());

        germlineSVs = getGermlineSVs();
        assertEquals(1, germlineSVs.size());
        assertTrue(germlineSVs.get(0).Reported);

        mLinx.clearClustersAndSVs();

        // intronic
        var = createSgl(mLinx.nextVarId(), CHR_1, 11500,POS_ORIENT);
        var.getSglMappings().add(new SglMapping(CHR_1, 13500, NEG_ORIENT, "", 1));

        mLinx.AllVariants.add(var);

        mLinx.preClusteringInit();
        mLinx.Analyser.clusterAndAnalyse();

        mDisruptionFinder.findReportableDisruptions(mLinx.AllVariants, mLinx.Analyser.getClusters());

        germlineSVs = getGermlineSVs();
        assertEquals(1, germlineSVs.size());
        assertTrue(germlineSVs.get(0).Reported);

        mLinx.clearClustersAndSVs();

        // same intronic is not disruptive
        var = createSgl(mLinx.nextVarId(), CHR_1, 11500,POS_ORIENT);
        var.getSglMappings().add(new SglMapping(CHR_1, 11800, NEG_ORIENT, "", 1));

        mLinx.AllVariants.add(var);

        mLinx.preClusteringInit();
        mLinx.Analyser.clusterAndAnalyse();

        mDisruptionFinder.findReportableDisruptions(mLinx.AllVariants, mLinx.Analyser.getClusters());

        germlineSVs = getGermlineSVs();
        assertTrue(germlineSVs.isEmpty());

        // evaluate on its own even if clustered with something else
        mLinx.clearClustersAndSVs();

        var = createSgl(mLinx.nextVarId(), CHR_1, 11500,POS_ORIENT);
        var.getSglMappings().add(new SglMapping(CHR_1, 13500, NEG_ORIENT, "", 1));

        SvVarData var2 = createSgl(mLinx.nextVarId(), CHR_1, 15000,POS_ORIENT);

        mLinx.AllVariants.add(var);
        mLinx.AllVariants.add(var2);

        mLinx.preClusteringInit();
        mLinx.Analyser.clusterAndAnalyse();

        assertEquals(1, mLinx.Analyser.getClusters().size());

        mDisruptionFinder.findReportableDisruptions(mLinx.AllVariants, mLinx.Analyser.getClusters());

        germlineSVs = getGermlineSVs();
        assertEquals(1, germlineSVs.size());
        assertTrue(germlineSVs.get(0).Reported);

        // ignore pseudogene
    }

    @Test
    public void testSglMappedDup()
    {
        // SGL with mapping to make a DUP around part of a gene
        SvVarData var = createSgl(mLinx.nextVarId(), CHR_1, 13500, POS_ORIENT);
        var.getSglMappings().add(new SglMapping(CHR_1, 11500, NEG_ORIENT, "", 1));

        mLinx.AllVariants.add(var);

        mLinx.preClusteringInit();
        mLinx.Analyser.clusterAndAnalyse();

        mDisruptionFinder.findReportableDisruptions(mLinx.AllVariants, mLinx.Analyser.getClusters());

        List<LinxGermlineSv> germlineSVs = getGermlineSVs();
        assertEquals(1, germlineSVs.size());
        assertTrue(germlineSVs.get(0).Reported);

        mLinx.clearClustersAndSVs();

        // cannot just be around the starting exon(s) of the gene
        var = createSgl(mLinx.nextVarId(), CHR_1, 13500, POS_ORIENT);
        var.getSglMappings().add(new SglMapping(CHR_1, 9000, NEG_ORIENT, "", 1));

        mLinx.AllVariants.add(var);

        mLinx.preClusteringInit();
        mLinx.Analyser.clusterAndAnalyse();

        mDisruptionFinder.findReportableDisruptions(mLinx.AllVariants, mLinx.Analyser.getClusters());

        germlineSVs = getGermlineSVs();
        assertTrue(germlineSVs.isEmpty());
    }


    private List<LinxGermlineSv> getGermlineSVs()
    {
        List<SvDisruptionData> standardDisruptions = Lists.newArrayList();
        List<LinxGermlineSv> germlineSVs = Lists.newArrayList();
        List<DriverCatalog> drivers = Lists.newArrayList();

        mDisruptionFinder.germlineDisruptions().populateGermlineSVs(standardDisruptions, germlineSVs, drivers);
        return germlineSVs;
    }

    private void addTestGeneData(EnsemblDataCache geneTransCache, final String chromosome, final String geneId, final String geneName)
    {
        byte strand = POS_STRAND;

        List<GeneData> geneList = Lists.newArrayList();
        geneList.add(createEnsemblGeneData(geneId, geneName, chromosome, strand, 10000, 120000));
        addGeneData(geneTransCache, chromosome, geneList);

        List<TranscriptData> transDataList = Lists.newArrayList();

        int transId = 1;

        int[] exonStarts = new int[] { 10000, 11000, 12000, 13000, 14000 };

        int codingStart = 11100;
        int codingEnd = 14100;
        TranscriptData transData = createTransExons(geneId, transId++, strand, exonStarts,  200, codingStart, codingEnd, true, "");
        transDataList.add(transData);

        addTransExonData(geneTransCache, geneId, transDataList);
    }
}
