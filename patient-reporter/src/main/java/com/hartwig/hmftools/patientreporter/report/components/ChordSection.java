package com.hartwig.hmftools.patientreporter.report.components;

import java.awt.Color;

import com.hartwig.hmftools.common.purple.purity.FittedPurityStatus;
import com.hartwig.hmftools.patientreporter.report.util.PatientReportFormat;

import org.jetbrains.annotations.NotNull;

import net.sf.dynamicreports.report.builder.component.ComponentBuilder;

public final class ChordSection {

    private static final int BUFFER = 3;
    private static final double START = 0;
    private static final double END = 1;

    @NotNull
    public static ComponentBuilder<?, ?> build(final double chordHrdProbability, @NotNull FittedPurityStatus fitStatus) {
        final int graphValue = computeGraphValue(chordHrdProbability);

//        final GradientBar gradient = ImmutableGradientBar.of(new Color(214, 234, 248), new Color(174, 214, 241), "Low", "High", graphValue);
        final GradientBar gradient = ImmutableGradientBar.of(new Color(239, 239, 239), new Color(171, 191, 171), "Low", "High", graphValue);

        final SliderSection sliderSection = ImmutableSliderSection.of("HR-Deficiency Probability",
                interpret(chordHrdProbability, fitStatus),
                description(),
                gradient);

        return sliderSection.build();
    }

    @NotNull
    private static String interpret(final double chordHrdProbability, @NotNull FittedPurityStatus fitStatus) {
        String hrdProbability = PatientReportFormat.formatPercent(chordHrdProbability);
        return PatientReportFormat.correctValueForFitStatus(fitStatus, hrdProbability);
    }

    private static int computeGraphValue(final double value) {
        final double scaledStart = scale(START);
        final double scaledEnd = scale(END);
        final double scaledIntervalLength = scaledEnd - scaledStart;
        final double scaleValue = Math.min(scaledEnd, Math.max(scaledStart, scale(value)));

        final double graphIntervalLength = 100 - 2 * BUFFER;
        return (int) Math.round((scaleValue - scaledStart) * graphIntervalLength / scaledIntervalLength + BUFFER);
    }

    private static double scale(double value) {
        return value;
    }

    @NotNull
    private static String description() {
        return "The HR-deficiency probability is determined by CHORD, a WGS signature-based classifier comparing the " +
                "signature of this sample with signatures found across samples with known BRCA1/BRCA2 inactivation.";
    }
}
