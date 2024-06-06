package com.hartwig.hmftools.orange;

import java.time.LocalDate;
import java.util.List;
import java.util.Map;
import java.util.Set;

import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import com.google.common.io.Resources;
import com.hartwig.hmftools.common.chord.ChordTestFactory;
import com.hartwig.hmftools.common.doid.DoidTestFactory;
import com.hartwig.hmftools.common.flagstat.FlagstatTestFactory;
import com.hartwig.hmftools.common.lilac.LilacTestFactory;
import com.hartwig.hmftools.common.linx.LinxTestFactory;
import com.hartwig.hmftools.common.metrics.WGSMetricsTestFactory;
import com.hartwig.hmftools.common.peach.PeachTestFactory;
import com.hartwig.hmftools.datamodel.hla.ImmutableLilacRecord;
import com.hartwig.hmftools.datamodel.hla.LilacAllele;
import com.hartwig.hmftools.datamodel.hla.LilacRecord;
import com.hartwig.hmftools.datamodel.immuno.ImmuneEscapeRecord;
import com.hartwig.hmftools.datamodel.isofox.AltSpliceJunctionContext;
import com.hartwig.hmftools.datamodel.isofox.AltSpliceJunctionType;
import com.hartwig.hmftools.datamodel.isofox.GeneExpression;
import com.hartwig.hmftools.datamodel.isofox.ImmutableIsofoxRecord;
import com.hartwig.hmftools.datamodel.isofox.IsofoxRecord;
import com.hartwig.hmftools.datamodel.isofox.IsofoxRnaStatistics;
import com.hartwig.hmftools.datamodel.isofox.NovelSpliceJunction;
import com.hartwig.hmftools.datamodel.isofox.RnaFusion;
import com.hartwig.hmftools.datamodel.isofox.StructuralVariantType;
import com.hartwig.hmftools.datamodel.linx.ImmutableLinxRecord;
import com.hartwig.hmftools.datamodel.linx.LinxFusion;
import com.hartwig.hmftools.datamodel.linx.LinxRecord;
import com.hartwig.hmftools.datamodel.orange.ExperimentType;
import com.hartwig.hmftools.datamodel.orange.ImmutableOrangePlots;
import com.hartwig.hmftools.datamodel.orange.ImmutableOrangeRecord;
import com.hartwig.hmftools.datamodel.orange.ImmutableOrangeSample;
import com.hartwig.hmftools.datamodel.orange.OrangePlots;
import com.hartwig.hmftools.datamodel.orange.OrangeRecord;
import com.hartwig.hmftools.datamodel.orange.OrangeRefGenomeVersion;
import com.hartwig.hmftools.datamodel.orange.OrangeSample;
import com.hartwig.hmftools.datamodel.peach.PeachGenotype;
import com.hartwig.hmftools.datamodel.purple.ImmutablePurpleRecord;
import com.hartwig.hmftools.datamodel.purple.PurpleRecord;
import com.hartwig.hmftools.datamodel.virus.ImmutableVirusInterpreterData;
import com.hartwig.hmftools.datamodel.virus.VirusBreakendQCStatus;
import com.hartwig.hmftools.datamodel.virus.VirusInterpretation;
import com.hartwig.hmftools.datamodel.virus.VirusInterpreterData;
import com.hartwig.hmftools.datamodel.virus.VirusInterpreterEntry;
import com.hartwig.hmftools.datamodel.virus.VirusLikelihoodType;
import com.hartwig.hmftools.orange.algo.cuppa.TestCuppaFactory;
import com.hartwig.hmftools.orange.algo.immuno.TestImmuneEscapeFactory;
import com.hartwig.hmftools.orange.algo.isofox.OrangeIsofoxTestFactory;
import com.hartwig.hmftools.orange.algo.linx.TestLinxInterpretationFactory;
import com.hartwig.hmftools.orange.algo.purple.TestPurpleInterpretationFactory;
import com.hartwig.hmftools.orange.algo.purple.TestPurpleVariantFactory;
import com.hartwig.hmftools.orange.algo.wildtype.TestWildTypeFactory;
import com.hartwig.hmftools.orange.conversion.LinxConversion;
import com.hartwig.hmftools.orange.conversion.OrangeConversion;

import org.apache.logging.log4j.util.Strings;
import org.jetbrains.annotations.NotNull;

public final class TestOrangeReportFactory
{
    private static final String TEST_SAMPLE = "TEST";
    private static final String DUMMY_IMAGE = Resources.getResource("test_images/white.png").getPath();

    @NotNull
    public static ImmutableOrangeRecord.Builder builder()
    {
        return ImmutableOrangeRecord.builder()
                .sampleId(TEST_SAMPLE)
                .samplingDate(LocalDate.of(2021, 11, 19))
                .experimentType(ExperimentType.TARGETED)
                .refGenomeVersion(OrangeRefGenomeVersion.V37)
                .tumorSample(createMinimalOrangeSample())
                .purple(TestPurpleInterpretationFactory.createMinimalTestPurpleData())
                .linx(TestLinxInterpretationFactory.createMinimalTestLinxData())
                .lilac(ImmutableLilacRecord.builder().qc(Strings.EMPTY).build())
                .immuneEscape(TestImmuneEscapeFactory.builder().build())
                .virusInterpreter(ImmutableVirusInterpreterData.builder().build())
                .chord(OrangeConversion.convert(ChordTestFactory.createMinimalTestChordAnalysis()))
                .cuppa(TestCuppaFactory.createMinimalCuppaData())
                .plots(createMinimalOrangePlots());
    }

    @NotNull
    public static OrangeRecord createMinimalTestReport()
    {
        return builder().build();
    }

    @NotNull
    public static OrangeRecord createProperTestReport()
    {
        return builder().experimentType(ExperimentType.WHOLE_GENOME)
                .addConfiguredPrimaryTumor(OrangeConversion.convert(DoidTestFactory.createDoidNode("1", "cancer type")))
                .platinumVersion("v5.35")
                .refSample(createMinimalOrangeSample())
                .germlineMVLHPerGene(createTestGermlineMVLHPerGene())
                .purple(createTestPurpleData())
                .linx(createTestLinxData())
                .addWildTypeGenes(TestWildTypeFactory.create("gene"))
                .isofox(createTestIsofoxData())
                .lilac(createTestLilacData())
                .immuneEscape(createTestImmuneEscapeRecord())
                .virusInterpreter(createTestVirusInterpreterData())
                .peach(createTestPeachData())
                .build();
    }

    @NotNull
    private static OrangeSample createMinimalOrangeSample()
    {
        return ImmutableOrangeSample.builder()
                .metrics(OrangeConversion.convert(WGSMetricsTestFactory.createMinimalTestWGSMetrics()))
                .flagstat(OrangeConversion.convert(FlagstatTestFactory.createMinimalTestFlagstat()))
                .build();
    }

    @NotNull
    private static OrangePlots createMinimalOrangePlots()
    {
        return ImmutableOrangePlots.builder()
                .sageTumorBQRPlot(DUMMY_IMAGE)
                .purpleInputPlot(DUMMY_IMAGE)
                .purpleFinalCircosPlot(DUMMY_IMAGE)
                .purpleClonalityPlot(DUMMY_IMAGE)
                .purpleCopyNumberPlot(DUMMY_IMAGE)
                .purpleVariantCopyNumberPlot(DUMMY_IMAGE)
                .purplePurityRangePlot(DUMMY_IMAGE)
                .purpleKataegisPlot(DUMMY_IMAGE)
                .build();
    }

    @NotNull
    private static Map<String, Double> createTestGermlineMVLHPerGene()
    {
        Map<String, Double> germlineMVLHPerGene = Maps.newHashMap();
        germlineMVLHPerGene.put("gene", 0.01);
        return germlineMVLHPerGene;
    }

    @NotNull
    private static PurpleRecord createTestPurpleData()
    {
        return ImmutablePurpleRecord.builder()
                .from(TestPurpleInterpretationFactory.createMinimalTestPurpleData())
                .addReportableSomaticVariants(TestPurpleVariantFactory.builder()
                        .gene("ARID1A")
                        .canonicalImpact(TestPurpleVariantFactory.impactBuilder()
                                .hgvsCodingImpact("c.1920+9571_1920+9596delAGTGAACCGTTGACTAGAGTTTGGTT")
                                .build())
                        .build())
                .addReportableSomaticVariants(TestPurpleVariantFactory.builder()
                        .gene("USH2A")
                        .canonicalImpact(TestPurpleVariantFactory.impactBuilder()
                                .hgvsCodingImpact("c.8558+420_8558+442delCCGATACGATGAAAGAAAAGAGC")
                                .build())
                        .build())
                .addReportableSomaticVariants(TestPurpleVariantFactory.builder()
                        .gene("USH2A")
                        .canonicalImpact(TestPurpleVariantFactory.impactBuilder().hgvsCodingImpact("c.11712-884A>T").build())
                        .addLocalPhaseSets(42256)
                        .build())
                .allGermlineVariants(Lists.newArrayList())
                .reportableGermlineVariants(Lists.newArrayList())
                .additionalSuspectGermlineVariants(Lists.newArrayList())
                .reportableGermlineFullLosses(Lists.newArrayList())
                .reportableGermlineLossOfHeterozygosities(Lists.newArrayList())
                .build();
    }

    @NotNull
    private static LinxRecord createTestLinxData()
    {
        LinxFusion fusion = LinxConversion.convert(LinxTestFactory.createMinimalTestFusion());
        return ImmutableLinxRecord.builder()
                .from(TestLinxInterpretationFactory.createMinimalTestLinxData())
                .addReportableSomaticFusions(fusion)
                .addReportableSomaticFusions(fusion)
                .addReportableSomaticFusions(fusion)
                .addReportableSomaticFusions(fusion)
                .addReportableSomaticFusions(fusion)
                .addReportableSomaticFusions(fusion)
                .addReportableSomaticFusions(fusion)
                .addReportableSomaticFusions(fusion)
                .allGermlineStructuralVariants(Lists.newArrayList())
                .allGermlineBreakends(Lists.newArrayList())
                .reportableGermlineBreakends(Lists.newArrayList())
                .build();
    }

    @NotNull
    private static LilacRecord createTestLilacData()
    {
        List<LilacAllele> alleles = Lists.newArrayList();
        alleles.add(OrangeConversion.convert(LilacTestFactory.alleleBuilder().allele("Allele 1").build(), true, true));
        alleles.add(OrangeConversion.convert(LilacTestFactory.alleleBuilder()
                .allele("Allele 2")
                .somaticInframeIndel(1D)
                .build(), true, true));

        return ImmutableLilacRecord.builder().qc("PASS").alleles(alleles).build();
    }

    @NotNull
    private static ImmuneEscapeRecord createTestImmuneEscapeRecord()
    {
        return TestImmuneEscapeFactory.builder().hasHlaEscape(true).build();
    }

    @NotNull
    private static IsofoxRecord createTestIsofoxData()
    {
        IsofoxRnaStatistics statistics =
                OrangeIsofoxTestFactory.rnaStatisticsBuilder().totalFragments(120000).duplicateFragments(60000).build();

        GeneExpression highExpression = OrangeIsofoxTestFactory.geneExpressionBuilder()
                .gene("MYC")
                .tpm(126.27)
                .medianTpmCancer(41D)
                .percentileCancer(0.91)
                .medianTpmCohort(37D)
                .percentileCohort(0.93)
                .build();

        GeneExpression lowExpression = OrangeIsofoxTestFactory.geneExpressionBuilder()
                .gene("CDKN2A")
                .tpm(5.34)
                .medianTpmCancer(18.32)
                .percentileCancer(0.04)
                .medianTpmCohort(16D)
                .percentileCohort(0.07)
                .build();

        RnaFusion novelKnownFusion = OrangeIsofoxTestFactory.rnaFusionBuilder()
                .geneStart("PTPRK")
                .geneEnd("RSPO3")
                .chromosomeStart("6")
                .positionStart(128841405)
                .chromosomeEnd("6")
                .positionEnd(127469792)
                .svType(StructuralVariantType.INV)
                .junctionTypeEnd("KNOWN")
                .junctionTypeEnd("KNOWN")
                .depthStart(73)
                .depthEnd(49)
                .splitFragments(8)
                .realignedFrags(0)
                .discordantFrags(1)
                .cohortFrequency(3)
                .build();

        RnaFusion novelPromiscuousFusion = OrangeIsofoxTestFactory.rnaFusionBuilder()
                .geneStart("NAP1L4")
                .geneEnd("BRAF")
                .chromosomeStart("11")
                .positionStart(2972480)
                .chromosomeEnd("7")
                .positionEnd(140487380)
                .svType(StructuralVariantType.BND)
                .junctionTypeStart("KNOWN")
                .junctionTypeEnd("KNOWN")
                .depthStart(9)
                .depthEnd(19)
                .splitFragments(5)
                .realignedFrags(2)
                .discordantFrags(3)
                .cohortFrequency(1)
                .build();

        NovelSpliceJunction novelSkippedExon = OrangeIsofoxTestFactory.novelSpliceJunctionBuilder()
                .chromosome("1")
                .gene("ALK")
                .junctionStart(50403003)
                .junctionEnd(60403003)
                .type(AltSpliceJunctionType.SKIPPED_EXONS)
                .depthStart(12)
                .depthEnd(14)
                .regionStart(AltSpliceJunctionContext.SPLICE_JUNC)
                .regionEnd(AltSpliceJunctionContext.SPLICE_JUNC)
                .fragmentCount(5)
                .cohortFrequency(3)
                .build();

        NovelSpliceJunction novelIntron = OrangeIsofoxTestFactory.novelSpliceJunctionBuilder()
                .chromosome("1")
                .gene("ALK")
                .junctionStart(50403003)
                .junctionEnd(60403003)
                .type(AltSpliceJunctionType.NOVEL_INTRON)
                .depthStart(24)
                .depthEnd(12)
                .regionStart(AltSpliceJunctionContext.EXONIC)
                .regionEnd(AltSpliceJunctionContext.EXONIC)
                .fragmentCount(43)
                .cohortFrequency(12)
                .build();

        return ImmutableIsofoxRecord.builder()
                .summary(statistics)
                .addAllAllGeneExpressions(Lists.newArrayList(highExpression, lowExpression))
                .addReportableHighExpression(highExpression)
                .addReportableLowExpression(lowExpression)
                .addAllAllFusions(Lists.newArrayList(novelKnownFusion, novelPromiscuousFusion))
                .addReportableNovelKnownFusions(novelKnownFusion)
                .addReportableNovelPromiscuousFusions(novelPromiscuousFusion)
                .addAllAllNovelSpliceJunctions(Lists.newArrayList(novelSkippedExon, novelIntron))
                .addReportableSkippedExons(novelSkippedExon)
                .addReportableNovelExonsIntrons(novelIntron)
                .build();
    }

    @NotNull
    private static VirusInterpreterData createTestVirusInterpreterData()
    {
        List<VirusInterpreterEntry> reportableViruses = Lists.newArrayList();

        reportableViruses.add(com.hartwig.hmftools.datamodel.virus.ImmutableVirusInterpreterEntry.builder()
                .name("virus A")
                .qcStatus(VirusBreakendQCStatus.NO_ABNORMALITIES)
                .integrations(3)
                .interpretation(VirusInterpretation.HPV)
                .percentageCovered(87D)
                .meanCoverage(42D)
                .expectedClonalCoverage(3D)
                .reported(true)
                .driverLikelihood(VirusLikelihoodType.UNKNOWN)
                .build());

        return ImmutableVirusInterpreterData.builder().reportableViruses(reportableViruses).build();
    }

    @NotNull
    private static Set<PeachGenotype> createTestPeachData()
    {
        return Set.of(OrangeConversion.convert(PeachTestFactory.builder().gene("DPYD").allele("allele").build()));
    }
}
