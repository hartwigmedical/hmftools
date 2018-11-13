package com.hartwig.hmftools.patientreporter.report.components;

import java.awt.Color;
import java.text.DecimalFormat;

import com.hartwig.hmftools.patientdb.data.Patient;
import com.hartwig.hmftools.patientreporter.report.util.PatientReportFormat;

import org.jetbrains.annotations.NotNull;

import net.sf.dynamicreports.report.builder.component.ComponentBuilder;

public final class MicrosatelliteSection {
    public static final double MSI_THRESHOLD = 4D;

    private static final int BUFFER = 3;
    private static final double START = 1E-2;
    private static final double END = 100;
    private static final Color COLOUR_BEGIN = new Color(239, 239, 239);
    private static final Color COLOUR_END = new Color(171, 191, 171);

    @NotNull
    public static ComponentBuilder<?, ?> build(double microsatelliteIndicator, boolean hasReliablePurityFit) {
        final int graphValue = computeGraphValue(microsatelliteIndicator);
        final int markerValue = computeGraphValue(MSI_THRESHOLD);

        final GradientBar gradient = !hasReliablePurityFit
                ? ImmutableGradientBar.ofOnlyMarker(COLOUR_BEGIN, COLOUR_END, "MSS", "MSI", markerValue)
                : ImmutableGradientBar.of(COLOUR_BEGIN, COLOUR_END, "MSS", "MSI", graphValue, markerValue);
        final SliderSection sliderSection = ImmutableSliderSection.of("Microsatellite Status",
                interpret(microsatelliteIndicator, hasReliablePurityFit),
                description(),
                gradient);
        return sliderSection.build();
    }

    @NotNull
    private static String interpret(double microsatelliteIndicator, boolean hasReliablePurityFit) {
        String formattedMs = new DecimalFormat("#.####").format(microsatelliteIndicator);
        String interpretedMs;
        
        if (microsatelliteIndicator > MSI_THRESHOLD) {
            interpretedMs = "Unstable (" + formattedMs + ")";
        } else {
            interpretedMs = "Stable (" + formattedMs + ")";
        }

        return PatientReportFormat.correctValueForFitReliability(interpretedMs, hasReliablePurityFit);
    }

    private static int computeGraphValue(double value) {
        final double scaledStart = scale(START);
        final double scaledEnd = scale(END);
        final double scaledValue = scale(value);
        final double scaledIntervalLength = scaledEnd - scaledStart;

        final int graphIntervalLength = 100 - 2 * BUFFER;
        return (int) Math.round((scaledValue - scaledStart) * graphIntervalLength / scaledIntervalLength + BUFFER);
    }

    private static double scale(double value) {
        return Math.log10(value);
    }

    @NotNull
    private static String description() {
        return "The microsatellite stability score represents the number of somatic inserts and deletes in (short) repeat sections "
                + "across the whole genome of the tumor per Mb. This metric can be considered as a good marker for instability "
                + "in microsatellite repeat regions. Tumors with a score greater than " + MSI_THRESHOLD + " are considered microsatellite "
                + "unstable (MSI).";
    }
}
