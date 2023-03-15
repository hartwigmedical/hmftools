package com.hartwig.hmftools.datamodel;

import com.hartwig.hmftools.datamodel.chord.ChordRecord;
import com.hartwig.hmftools.datamodel.chord.ChordStatus;
import com.hartwig.hmftools.datamodel.chord.TestChordFactory;
import com.hartwig.hmftools.datamodel.cuppa.CuppaData;
import com.hartwig.hmftools.datamodel.cuppa.ImmutableCuppaData;
import com.hartwig.hmftools.datamodel.cuppa.TestCuppaFactory;
import com.hartwig.hmftools.datamodel.hla.ImmutableLilacRecord;
import com.hartwig.hmftools.datamodel.hla.LilacRecord;
import com.hartwig.hmftools.datamodel.lilac.TestLilacFactory;
import com.hartwig.hmftools.datamodel.linx.*;
import com.hartwig.hmftools.datamodel.orange.*;
import com.hartwig.hmftools.datamodel.peach.PeachGenotype;
import com.hartwig.hmftools.datamodel.peach.TestPeachFactory;
import com.hartwig.hmftools.datamodel.purple.*;
import com.hartwig.hmftools.datamodel.sv.LinxBreakendType;
import com.hartwig.hmftools.datamodel.virus.*;
import org.apache.logging.log4j.util.Strings;
import org.jetbrains.annotations.NotNull;

import java.time.LocalDate;
import java.util.Set;

public final class TestOrangeFactory {

    private static final String EMPTY_CIRCOS_PLOT = Thread.currentThread().getContextClassLoader().getResource("plot/empty.circos.png").getPath();

    private TestOrangeFactory() {
    }

    @NotNull
    public static OrangeRecord createMinimalTestOrangeRecord() {
        return ImmutableOrangeRecord.builder()
                .sampleId("TEST")
                .experimentDate(LocalDate.of(2022, 1, 20))
                .refGenomeVersion(OrangeRefGenomeVersion.V37)
                .purple(createMinimalTestPurpleRecord())
                .linx(ImmutableLinxRecord.builder().build())
                .peach(Set.of())
                .cuppa(createTestCuppaRecord())
                .virusInterpreter(ImmutableVirusInterpreterData.builder().build())
                .lilac(createMinimalTestLilacRecord())
                .chord(createTestChordRecord())
                .plots(createMinimalTestOrangePlots())
                .build();
    }

    @NotNull
    private static PurpleRecord createMinimalTestPurpleRecord() {
        return ImmutablePurpleRecord.builder()
                .fit(TestPurpleFactory.fitBuilder().build())
                .characteristics(TestPurpleFactory.characteristicsBuilder().build())
                .build();
    }

    @NotNull
    private static LilacRecord createMinimalTestLilacRecord() {
        return ImmutableLilacRecord.builder().qc(Strings.EMPTY).build();
    }

    @NotNull
    private static OrangePlots createMinimalTestOrangePlots() {
        return ImmutableOrangePlots.builder()
                .sageTumorBQRPlot(Strings.EMPTY)
                .purpleInputPlot(Strings.EMPTY)
                .purpleClonalityPlot(Strings.EMPTY)
                .purpleCopyNumberPlot(Strings.EMPTY)
                .purpleVariantCopyNumberPlot(Strings.EMPTY)
                .purplePurityRangePlot(Strings.EMPTY)
                .purpleFinalCircosPlot(EMPTY_CIRCOS_PLOT)
                .build();
    }

    @NotNull
    public static OrangeRecord createProperTestOrangeRecord() {
        return ImmutableOrangeRecord.builder()
                .from(createMinimalTestOrangeRecord())
                .purple(createTestPurpleRecord())
                .linx(createTestLinxRecord())
                .peach(createTestPeachRecord())
                .cuppa(createTestCuppaRecord())
                .virusInterpreter(createTestVirusInterpreterRecord())
                .lilac(createTestLilacRecord())
                .chord(createTestChordRecord())
                .build();
    }

    @NotNull
    private static PurpleRecord createTestPurpleRecord() {
        PurpleVariant variant = TestPurpleFactory.variantBuilder()
                .reported(true)
                .gene("BRAF")
                .adjustedCopyNumber(6.0)
                .variantCopyNumber(4.1)
                .hotspot(Hotspot.HOTSPOT)
                .subclonalLikelihood(0.02)
                .biallelic(false)
                .canonicalImpact(TestPurpleFactory.transcriptImpactBuilder()
                        .hgvsCodingImpact("c.something")
                        .hgvsProteinImpact("p.Val600Glu")
                        .spliceRegion(false)
                        .addEffects(PurpleVariantEffect.MISSENSE)
                        .codingEffect(PurpleCodingEffect.MISSENSE)
                        .build())
                .build();

        PurpleCopyNumber copyNumber =
                TestPurpleFactory.copyNumberBuilder().chromosome("1").start(10).end(20).averageTumorCopyNumber(2.1).build();

        PurpleGainLoss gain = TestPurpleFactory.gainLossBuilder().gene("MYC").interpretation(CopyNumberInterpretation.FULL_GAIN)
                .minCopies(38)
                .maxCopies(40)
                .build();

        PurpleGainLoss loss = TestPurpleFactory.gainLossBuilder().gene("PTEN").interpretation(CopyNumberInterpretation.FULL_LOSS)
                .minCopies(0)
                .maxCopies(0)
                .build();

        return ImmutablePurpleRecord.builder()
                .from(createMinimalTestPurpleRecord())
                .fit(createTestPurpleFit())
                .characteristics(createTestPurpleCharacteristics())
                .addSomaticDrivers(TestPurpleFactory.driverBuilder()
                        .gene(variant.gene())
                        .driver(PurpleDriverType.MUTATION)
                        .driverLikelihood(1D)
                        .build())
                .addSomaticDrivers(TestPurpleFactory.driverBuilder()
                        .gene(gain.gene())
                        .driver(PurpleDriverType.AMP)
                        .driverLikelihood(1D)
                        .build())
                .addSomaticDrivers(TestPurpleFactory.driverBuilder()
                        .gene(loss.gene())
                        .driver(PurpleDriverType.DEL)
                        .driverLikelihood(1D)
                        .build())
                .addAllSomaticVariants(variant)
                .addReportableSomaticVariants(variant)
                .addAllSomaticCopyNumbers(copyNumber)
                .addAllSomaticGainsLosses(gain, loss)
                .addReportableSomaticGainsLosses(gain, loss)
                .build();
    }

    @NotNull
    private static PurpleFit createTestPurpleFit() {
        return TestPurpleFactory.fitBuilder()
                .hasSufficientQuality(true)
                .containsTumorCells(false)
                .purity(0.12)
                .ploidy(3.1)
                .build();
    }

    @NotNull
    private static PurpleCharacteristics createTestPurpleCharacteristics() {
        return TestPurpleFactory.characteristicsBuilder()
                .microsatelliteIndelsPerMb(0.1)
                .microsatelliteStatus(PurpleMicrosatelliteStatus.MSS)
                .tumorMutationalBurdenPerMb(13D)
                .tumorMutationalBurdenStatus(PurpleTumorMutationalStatus.HIGH)
                .tumorMutationalLoad(189)
                .tumorMutationalLoadStatus(PurpleTumorMutationalStatus.HIGH)
                .build();
    }

    @NotNull
    private static LinxRecord createTestLinxRecord() {
        LinxBreakend breakend1 = TestLinxFactory.breakendBuilder()
                .reportedDisruption(true)
                .svId(1)
                .gene("RB1")
                .type(LinxBreakendType.DEL)
                .junctionCopyNumber(0.8)
                .undisruptedCopyNumber(2.1)
                .build();

        LinxBreakend breakend2 = TestLinxFactory.breakendBuilder()
                .reportedDisruption(true)
                .svId(1)
                .gene("PTEN")
                .type(LinxBreakendType.DEL)
                .junctionCopyNumber(1D)
                .undisruptedCopyNumber(1D)
                .build();

        LinxFusion fusion = TestLinxFactory.fusionBuilder()
                .reported(true)
                .reportedType(LinxFusionType.KNOWN_PAIR)
                .geneStart("EML4")
                .fusedExonUp(2)
                .geneEnd("ALK")
                .fusedExonDown(4)
                .likelihood(FusionLikelihoodType.HIGH)
                .build();

        return ImmutableLinxRecord.builder()
                .addAllSomaticStructuralVariants(TestLinxFactory.structuralVariantBuilder().svId(1).clusterId(1).build())
                .addSomaticHomozygousDisruptions(TestLinxFactory.homozygousDisruptionBuilder().gene("TP53").build())
                .addAllSomaticBreakends(breakend1, breakend2)
                .addReportableSomaticBreakends(breakend1, breakend2)
                .addAllSomaticFusions(fusion)
                .addReportableSomaticFusions(fusion)
                .build();
    }

    @NotNull
    private static Set<PeachGenotype> createTestPeachRecord() {
        return Set.of(TestPeachFactory.builder()
                .gene("DPYD")
                .haplotype("1* HOM")
                .function("Normal function")
                .linkedDrugs("5-Fluorouracil")
                .urlPrescriptionInfo("https://www.pharmgkb.org/guidelineAnnotation/PA166104939")
                .panelVersion("peach_prod_v1.3")
                .repoVersion("1.7")
                .build());
    }

    @NotNull
    private static CuppaData createTestCuppaRecord() {
        return ImmutableCuppaData.builder()
                .addPredictions(TestCuppaFactory.builder().cancerType("Melanoma").likelihood(0.996).build())
                .simpleDups32To200B(0)
                .maxComplexSize(0)
                .telomericSGLs(0)
                .lineCount(0)
                .build();
    }

    @NotNull
    private static VirusInterpreterData createTestVirusInterpreterRecord() {
        AnnotatedVirus virus = TestVirusInterpreterFactory.builder()
                .reported(true)
                .name("Human papillomavirus type 16")
                .qcStatus(VirusBreakendQCStatus.NO_ABNORMALITIES)
                .interpretation(VirusInterpretation.HPV)
                .integrations(3)
                .virusDriverLikelihoodType(VirusLikelihoodType.HIGH)
                .meanCoverage(0)
                .build();

        return ImmutableVirusInterpreterData.builder().addAllViruses(virus).addReportableViruses(virus).build();
    }

    @NotNull
    private static LilacRecord createTestLilacRecord() {
        return ImmutableLilacRecord.builder()
                .qc("PASS")
                .addAlleles(TestLilacFactory.builder().allele("A*01:01").tumorCopyNumber(1.2).build())
                .build();
    }

    @NotNull
    private static ChordRecord createTestChordRecord() {
        return TestChordFactory.builder()
                .brca1Value(0)
                .brca2Value(0)
                .hrdType(Strings.EMPTY)
                .hrStatus(ChordStatus.HR_PROFICIENT).build();
    }
}
