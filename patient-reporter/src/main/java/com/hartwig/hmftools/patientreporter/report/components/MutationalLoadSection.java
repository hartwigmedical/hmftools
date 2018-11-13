package com.hartwig.hmftools.patientreporter.report.components;

import java.awt.Color;
import java.text.DecimalFormat;

import com.hartwig.hmftools.patientreporter.report.util.PatientReportFormat;

import org.jetbrains.annotations.NotNull;

import net.sf.dynamicreports.report.builder.component.ComponentBuilder;

public final class MutationalLoadSection {
    public static final int ML_THRESHOLD = 140;

    private static final int BUFFER = 3;
    private static final int START = 1;
    private static final int END = 1000;

    @NotNull
    public static ComponentBuilder<?, ?> build(int mutationalLoad, boolean hasReliablePurityFit) {

        final int graphValue = computeGraphValue(mutationalLoad);
        final int markerValue = computeGraphValue(ML_THRESHOLD);
        final GradientBar gradient = !hasReliablePurityFit ?
                ImmutableGradientBar.ofOnlyMarker(new Color(239, 229, 203), new Color(159, 163, 193), "Low", "High", markerValue) :
                ImmutableGradientBar.of(new Color(239, 229, 203), new Color(159, 163, 193), "Low", "High", graphValue, markerValue);


        final SliderSection sliderSection = ImmutableSliderSection.of("Tumor Mutational Load",
                interpret(mutationalLoad, hasReliablePurityFit),
                description(),
                gradient);
        return sliderSection.build();
    }

    @NotNull
    private static String interpret(int mutationalLoad, boolean hasReliablePurityFit) {
        String interpretedML;

        if (mutationalLoad > ML_THRESHOLD) {
            interpretedML = "High (" + mutationalLoad + ")";
        } else {
            interpretedML = "Low (" + mutationalLoad + ")";
        }

        return PatientReportFormat.correctValueForFitReliability(interpretedML, hasReliablePurityFit);

    }

    private static int computeGraphValue(double value) {
        final double scaledStart = scale(START);
        final double scaledEnd = scale(END);
        final double scaledIntervalLength = scaledEnd - scaledStart;
        final double scaleValue = Math.min(scaledEnd, Math.max(scaledStart, scale(value)));

        final double graphIntervalLength = 100 - 2 * BUFFER;
        return (int) Math.round((scaleValue - scaledStart) * graphIntervalLength / scaledIntervalLength + BUFFER);
    }

    private static double scale(double value) {
        return Math.log10(value);
    }

    @NotNull
    private static String description() {
        return "The tumor mutational load represents the total number of somatic missense variants across"
                + " the whole genome of the tumor. Patients with a mutational load over " + ML_THRESHOLD
                + " could be eligible for immunotherapy within the DRUP study. ";
    }
}
