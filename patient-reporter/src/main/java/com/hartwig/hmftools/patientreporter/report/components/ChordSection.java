package com.hartwig.hmftools.patientreporter.report.components;

import java.awt.Color;
import java.text.DecimalFormat;

import com.hartwig.hmftools.patientreporter.report.util.PatientReportFormat;

import org.jetbrains.annotations.NotNull;

import net.sf.dynamicreports.report.builder.component.ComponentBuilder;

public final class ChordSection {

    private static final int BUFFER = 3;
    private static final double START = 0;
    private static final double END = 1;

    @NotNull
    public static ComponentBuilder<?, ?> build(double chordHrdScore, boolean hasReliablePurityFit) {
        final int graphValue = computeGraphValue(chordHrdScore);

        final GradientBar gradient = ImmutableGradientBar.of(new Color(214, 234, 248), new Color(174, 214, 241), "Low", "High", graphValue);

        final SliderSection sliderSection =
                ImmutableSliderSection.of("HR-Deficiency Score", interpret(chordHrdScore, hasReliablePurityFit), description(), gradient);

        return sliderSection.build();
    }

    @NotNull
    private static String interpret(double chordHrdScore, boolean hasReliablePurityFit) {
        String hrdScoreString = new DecimalFormat("#.##").format(chordHrdScore);
        return PatientReportFormat.correctValueForFitReliability(hrdScoreString, hasReliablePurityFit);
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
        return "The HR-deficiency score is determined by CHORD, a WGS signature-based classifier comparing the "
                + "signature of this sample with signatures found across samples with known BRCA1/BRCA2 inactivation.";
    }
}
