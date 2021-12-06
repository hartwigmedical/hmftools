package com.hartwig.hmftools.orange;

import java.time.LocalDate;
import java.util.List;

import com.google.common.collect.Lists;
import com.google.common.io.Resources;
import com.hartwig.hmftools.common.chord.ChordTestFactory;
import com.hartwig.hmftools.common.cuppa.CuppaTestFactory;
import com.hartwig.hmftools.common.flagstat.FlagstatTestFactory;
import com.hartwig.hmftools.common.linx.ImmutableLinxData;
import com.hartwig.hmftools.common.linx.LinxData;
import com.hartwig.hmftools.common.linx.LinxTestFactory;
import com.hartwig.hmftools.common.metrics.WGSMetricsTestFactory;
import com.hartwig.hmftools.common.protect.ImmutableProtectEvidence;
import com.hartwig.hmftools.common.protect.ProtectEvidence;
import com.hartwig.hmftools.common.protect.ProtectTestFactory;
import com.hartwig.hmftools.common.purple.ImmutablePurpleData;
import com.hartwig.hmftools.common.purple.PurpleData;
import com.hartwig.hmftools.common.purple.PurpleTestFactory;
import com.hartwig.hmftools.common.sv.linx.ImmutableLinxFusion;
import com.hartwig.hmftools.common.sv.linx.LinxFusion;
import com.hartwig.hmftools.common.variant.ImmutableReportableVariant;
import com.hartwig.hmftools.common.variant.VariantTestFactory;
import com.hartwig.hmftools.common.virus.ImmutableVirusInterpreterData;
import com.hartwig.hmftools.orange.algo.ImmutableOrangePlots;
import com.hartwig.hmftools.orange.algo.ImmutableOrangeReport;
import com.hartwig.hmftools.orange.algo.ImmutableOrangeSample;
import com.hartwig.hmftools.orange.algo.OrangePlots;
import com.hartwig.hmftools.orange.algo.OrangeReport;
import com.hartwig.hmftools.orange.algo.OrangeSample;

import org.jetbrains.annotations.NotNull;

public final class OrangeReportTestFactory {

    private static final String TEST_SAMPLE = "TEST";
    private static final String DUMMY_IMAGE = Resources.getResource("test_images/white.png").getPath();

    private OrangeReportTestFactory() {
    }

    @NotNull
    public static OrangeReport createMinimalTestReport() {
        return ImmutableOrangeReport.builder()
                .sampleId(TEST_SAMPLE)
                .reportDate(LocalDate.of(2021, 11, 19))
                .refSample(createMinimalOrangeSample())
                .tumorSample(createMinimalOrangeSample())
                .purple(PurpleTestFactory.createMinimalTestPurpleData())
                .linx(ImmutableLinxData.builder().build())
                .virusInterpreter(ImmutableVirusInterpreterData.builder().build())
                .chord(ChordTestFactory.createMinimalTestChordAnalysis())
                .cuppa(CuppaTestFactory.createMinimalCuppaData())
                .plots(createMinimalOrangePlots())
                .build();
    }

    @NotNull
    public static OrangeReport createProperTestReport() {
        return ImmutableOrangeReport.builder()
                .from(createMinimalTestReport())
                .purple(createTestPurpleData())
                .linx(createTestLinxData())
                .protect(createTestProtectData())
                .build();
    }

    @NotNull
    private static OrangeSample createMinimalOrangeSample() {
        return ImmutableOrangeSample.builder()
                .metrics(WGSMetricsTestFactory.createMinimalTestWGSMetrics())
                .flagstat(FlagstatTestFactory.createMinimalTestFlagstat())
                .build();
    }

    @NotNull
    private static OrangePlots createMinimalOrangePlots() {
        return ImmutableOrangePlots.builder()
                .sageReferenceBQRPlot(DUMMY_IMAGE)
                .sageTumorBQRPlot(DUMMY_IMAGE)
                .purpleInputPlot(DUMMY_IMAGE)
                .purpleFinalCircosPlot(DUMMY_IMAGE)
                .purpleClonalityPlot(DUMMY_IMAGE)
                .purpleCopyNumberPlot(DUMMY_IMAGE)
                .purpleVariantCopyNumberPlot(DUMMY_IMAGE)
                .purplePurityRangePlot(DUMMY_IMAGE)
                .cuppaSummaryPlot(DUMMY_IMAGE)
                .build();
    }

    @NotNull
    private static PurpleData createTestPurpleData() {
        return ImmutablePurpleData.builder()
                .from(PurpleTestFactory.createMinimalTestPurpleData())
                .addReportableSomaticVariants(ImmutableReportableVariant.builder()
                        .from(VariantTestFactory.createTestReportableVariant())
                        .gene("ARID1A")
                        .canonicalHgvsCodingImpact("c.1920+9571_1920+9596delAGTGAACCGTTGACTAGAGTTTGGTT")
                        .build())
                .addReportableSomaticVariants(ImmutableReportableVariant.builder()
                        .from(VariantTestFactory.createTestReportableVariant())
                        .gene("USH2A")
                        .canonicalHgvsCodingImpact("c.8558+420_8558+442delCCGATACGATGAAAGAAAAGAGC")
                        .build())
                .addReportableSomaticVariants(ImmutableReportableVariant.builder()
                        .from(VariantTestFactory.createTestReportableVariant())
                        .gene("USH2A")
                        .canonicalHgvsCodingImpact("c.11712-884A>T")
                        .localPhaseSet(42256)
                        .build())
                .build();
    }

    @NotNull
    private static LinxData createTestLinxData() {
        LinxFusion fusion = ImmutableLinxFusion.builder().from(LinxTestFactory.createMinimalTestFusion()).build();
        return ImmutableLinxData.builder()
                .addReportableFusions(fusion)
                .addReportableFusions(fusion)
                .addReportableFusions(fusion)
                .addReportableFusions(fusion)
                .addReportableFusions(fusion)
                .addReportableFusions(fusion)
                .addReportableFusions(fusion)
                .addReportableFusions(fusion)
                .build();
    }

    @NotNull
    private static List<ProtectEvidence> createTestProtectData() {
        List<ProtectEvidence> evidences = Lists.newArrayList();

        evidences.add(ImmutableProtectEvidence.builder()
                .from(ProtectTestFactory.createTestProtectEvidence())
                .genomicEvent("USH2A c.8558+420_8558+442delCCGATACGATGAAAGAAAAGAGC")
                .build());

        return evidences;
    }
}
