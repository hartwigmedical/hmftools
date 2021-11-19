package com.hartwig.hmftools.orange;

import java.time.LocalDate;

import com.google.common.io.Resources;
import com.hartwig.hmftools.common.chord.ChordTestFactory;
import com.hartwig.hmftools.common.cuppa.CuppaTestFactory;
import com.hartwig.hmftools.common.flagstat.FlagstatTestFactory;
import com.hartwig.hmftools.common.linx.ImmutableLinxData;
import com.hartwig.hmftools.common.metrics.WGSMetricsTestFactory;
import com.hartwig.hmftools.common.purple.ImmutablePurpleData;
import com.hartwig.hmftools.common.purple.PurpleData;
import com.hartwig.hmftools.common.purple.PurpleTestFactory;
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
        return ImmutableOrangeReport.builder().from(createMinimalTestReport()).purple(createTestPurpleData()).build();
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
                        .gene("USH2A")
                        .canonicalHgvsCodingImpact("c.8558+420_8558+442delCCGATACGATGAAAGAAAAGAGC")
                        .localPhaseSet(42256)
                        .build())
                .build();

    }
}
