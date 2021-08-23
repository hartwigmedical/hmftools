package com.hartwig.hmftools.linx.utils;

import static com.hartwig.hmftools.common.drivercatalog.DriverCategory.TSG;
import static com.hartwig.hmftools.common.test.GeneTestUtils.addGeneData;
import static com.hartwig.hmftools.common.test.GeneTestUtils.addTransExonData;
import static com.hartwig.hmftools.common.test.GeneTestUtils.createEnsemblGeneData;
import static com.hartwig.hmftools.common.test.GeneTestUtils.createTransExons;
import static com.hartwig.hmftools.common.test.GeneTestUtils.generateExonStarts;
import static com.hartwig.hmftools.common.gene.TranscriptProteinData.BIOTYPE_PROTEIN_CODING;
import static com.hartwig.hmftools.common.gene.CodingBaseData.PHASE_NONE;
import static com.hartwig.hmftools.common.fusion.FusionCommon.NEG_STRAND;
import static com.hartwig.hmftools.common.fusion.FusionCommon.POS_STRAND;

import java.util.List;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.drivercatalog.panel.DriverGene;
import com.hartwig.hmftools.common.drivercatalog.panel.DriverGeneGermlineReporting;
import com.hartwig.hmftools.common.drivercatalog.panel.ImmutableDriverGene;
import com.hartwig.hmftools.common.ensemblcache.EnsemblDataCache;
import com.hartwig.hmftools.common.gene.GeneData;
import com.hartwig.hmftools.common.gene.TranscriptData;
import com.hartwig.hmftools.common.fusion.BreakendGeneData;
import com.hartwig.hmftools.common.fusion.BreakendTransData;

public class GeneTestUtils
{
    public static final String CHR_1 = "1";

    public static final String GENE_NAME_1 = "GENE1"; // +ve strand
    public static final String GENE_ID_1 = "ENSG0001";
    public static final int GENE_START_1 = 1000;

    public static final String GENE_NAME_2 = "GENE2"; // +ve strand
    public static final String GENE_ID_2 = "ENSG0002";
    public static final int GENE_START_2 = 10000;

    public static final String GENE_NAME_3 = "GENE3"; // -ve strand
    public static final String GENE_ID_3 = "ENSG0003";
    public static final int GENE_START_3 = 20000;

    public static final String CHR_2 = "2";

    public static final String GENE_NAME_4 = "GENE4"; // +ve strand
    public static final String GENE_ID_4 = "ENSG0004";
    public static final int GENE_START_4 = 1000;

    public static final String GENE_NAME_5 = "GENE5"; // -ve strand
    public static final String GENE_ID_5 = "ENSG0005";
    public static final int GENE_START_5 = 10000;

    public static final String GENE_NAME_6 = "GENE6"; // -ve strand
    public static final String GENE_ID_6 = "ENSG0006";
    public static final int GENE_START_6 = 10400; // will overlap with previous gene and share some exons

    public static final int EXON_LENGTH = 100;

    public static void addTestGenes(EnsemblDataCache geneTransCache)
    {
        List<GeneData> geneList = Lists.newArrayList();
        geneList.add(createEnsemblGeneData(GENE_ID_1, GENE_NAME_1, CHR_1, POS_STRAND, GENE_START_1, GENE_START_1 + 1000));
        geneList.add(createEnsemblGeneData(GENE_ID_2, GENE_NAME_2, CHR_1, POS_STRAND, GENE_START_2, GENE_START_2 + 1000));
        geneList.add(createEnsemblGeneData(GENE_ID_3, GENE_NAME_3, CHR_1, NEG_STRAND, GENE_START_3, GENE_START_3 + 1000));
        addGeneData(geneTransCache, CHR_1, geneList);

        geneList = Lists.newArrayList();
        geneList.add(createEnsemblGeneData(GENE_ID_4, GENE_NAME_4, CHR_2, POS_STRAND, GENE_START_4, GENE_START_4 + 1000));
        geneList.add(createEnsemblGeneData(GENE_ID_5, GENE_NAME_5, CHR_2, NEG_STRAND, GENE_START_5, GENE_START_5 + 1500));
        geneList.add(createEnsemblGeneData(GENE_ID_6, GENE_NAME_6, CHR_2, NEG_STRAND, GENE_START_6, GENE_START_6 + 1500));
        addGeneData(geneTransCache, CHR_2, geneList);
    }

    public static final int TRANS_1 = 1;
    public static final int TRANS_2 = 2;
    public static final int TRANS_3 = 3;
    public static final int TRANS_4 = 4;
    public static final int TRANS_5 = 5;
    public static final int TRANS_6 = 5;

    public static void addTestTranscripts(EnsemblDataCache geneTransCache)
    {
        Integer codingStart = null;
        Integer codingEnd = null;
        boolean canonical = true;

        List<TranscriptData> transDataList = Lists.newArrayList();

        // exons will be gene start + 0-100, 200-300, 400-500

        TranscriptData transData = createTransExons(
                GENE_ID_1, TRANS_1, POS_STRAND, generateExonStarts(GENE_START_1, 7, EXON_LENGTH, 100),
                EXON_LENGTH, GENE_START_1 + 250, GENE_START_1 + 1450, canonical, BIOTYPE_PROTEIN_CODING);

        transDataList.add(transData);

        addTransExonData(geneTransCache, GENE_ID_1, transDataList);

        transDataList = Lists.newArrayList();

        codingStart = GENE_START_2 + 250;
        codingEnd = GENE_START_2 + 850;
        transData = createTransExons(
                GENE_ID_2, TRANS_2, POS_STRAND, generateExonStarts(GENE_START_2, 5, EXON_LENGTH, 100),
                EXON_LENGTH, GENE_START_2 + 250, GENE_START_2 + 850, canonical, BIOTYPE_PROTEIN_CODING);

        transDataList.add(transData);

        addTransExonData(geneTransCache, GENE_ID_2, transDataList);

        transDataList = Lists.newArrayList();

        transData = createTransExons(
                GENE_ID_3, TRANS_3, NEG_STRAND, generateExonStarts(GENE_START_3, 3, EXON_LENGTH, 100),
                EXON_LENGTH, GENE_START_3 + 250, GENE_START_3 + 850, canonical, BIOTYPE_PROTEIN_CODING);

        transDataList.add(transData);

        addTransExonData(geneTransCache, GENE_ID_3, transDataList);

        transDataList = Lists.newArrayList();

        transData = createTransExons(
                GENE_ID_4, TRANS_4, POS_STRAND, generateExonStarts(GENE_START_4, 3, EXON_LENGTH, 100),
                EXON_LENGTH, GENE_START_4 + 250, GENE_START_4 + 850, canonical, BIOTYPE_PROTEIN_CODING);

        transDataList.add(transData);

        addTransExonData(geneTransCache, GENE_ID_4, transDataList);

        transDataList = Lists.newArrayList();

        transData = createTransExons(
                GENE_ID_5, TRANS_5, NEG_STRAND, generateExonStarts(GENE_START_5, 5, EXON_LENGTH, 100),
                EXON_LENGTH, GENE_START_5 + 250, GENE_START_5 + 850, canonical, BIOTYPE_PROTEIN_CODING);

        transDataList.add(transData);

        addTransExonData(geneTransCache, GENE_ID_5, transDataList);

        transDataList = Lists.newArrayList();

        transData = createTransExons(
                GENE_ID_6, TRANS_6, NEG_STRAND, generateExonStarts(GENE_START_6, 6, EXON_LENGTH, 100),
                EXON_LENGTH, GENE_START_6 + 250, GENE_START_6 + 850, canonical, BIOTYPE_PROTEIN_CODING);

        transDataList.add(transData);

        addTransExonData(geneTransCache, GENE_ID_6, transDataList);
    }

    public static DriverGene createDriverGene(final String name)
    {
        return ImmutableDriverGene.builder()
                .gene(name)
                .reportMissenseAndInframe(false)
                .reportNonsenseAndFrameshift(false)
                .reportSplice(false)
                .reportDeletion(false)
                .reportDisruption(true)
                .reportAmplification(false)
                .reportSomaticHotspot(false)
                .reportGermlineVariant(DriverGeneGermlineReporting.NONE)
                .reportGermlineHotspot(DriverGeneGermlineReporting.NONE)
                .likelihoodType(TSG)
                .reportGermlineDisruption(false)
                .build();
    }
    public static BreakendTransData createTranscript(
            final BreakendGeneData gene, int transId, boolean isCanonical,
            int transStart, int transEnd, Integer codingStart, Integer codingEnd, String bioType,
            final int exonUpstream, final int exonDownstream, int phase, int codingBases, int totalCodingBases)
    {
        return createTranscript(gene, transId, isCanonical, transStart, transEnd, codingStart, codingEnd, bioType,
                exonUpstream, exonDownstream, phase, PHASE_NONE, codingBases, totalCodingBases);
    }

    public static BreakendTransData createTranscript(
            final BreakendGeneData gene, int transId, boolean isCanonical,
            int transStart, int transEnd, Integer codingStart, Integer codingEnd, String bioType,
            final int exonUpstream, final int exonDownstream, int phase, int exonicPhase, int codingBases, int totalCodingBases)
    {
        final String transName = String.format("TRANS%03d", transId);
        TranscriptData transData = new TranscriptData(
                transId, transName, gene.StableId, isCanonical, gene.Strand, transStart, transEnd, codingStart, codingEnd, bioType);

        return new BreakendTransData(gene, transData, exonUpstream, exonDownstream, phase, exonicPhase, codingBases, totalCodingBases);
    }

}
