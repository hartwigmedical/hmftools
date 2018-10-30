package com.hartwig.hmftools.patientreporter.report.components;

import java.awt.Color;
import java.text.DecimalFormat;

import com.hartwig.hmftools.patientreporter.report.util.PatientReportFormat;

import org.jetbrains.annotations.NotNull;

import net.sf.dynamicreports.report.builder.component.ComponentBuilder;

public final class MicrosatelliteSection {
    private static final double MSI_THRESHOLD = 4D;

    private static final int BUFFER = 3;
    private static final double START = 1E-2;
    private static final double END = 100;

    @NotNull
    public static ComponentBuilder<?, ?> build(double microsatelliteIndicator, boolean hasReliablePurityFit) {
        final int graphValue = computeGraphValue(microsatelliteIndicator);
        final int markerValue = computeGraphValue(MSI_THRESHOLD);

        final GradientBar gradient =
                ImmutableGradientBar.of(new Color(239, 239, 239), new Color(171, 191, 171), "MSS", "MSI", graphValue, markerValue);
        final SliderSection sliderSection = ImmutableSliderSection.of("Microsatellite Status",
                interpret(microsatelliteIndicator, hasReliablePurityFit),
                description(),
                gradient);
        return sliderSection.build();
    }

    @NotNull
    public static String interpret(double microsatelliteIndicator, boolean hasReliablePurityFit) {
        final String formattedMicrosatelliteIndicator =
                PatientReportFormat.correctValueForFitReliability(new DecimalFormat("#.####").format(microsatelliteIndicator),
                        hasReliablePurityFit);
        if (microsatelliteIndicator > MSI_THRESHOLD) {
            return "Unstable (" + formattedMicrosatelliteIndicator + ")";
        } else {
            return "Stable (" + formattedMicrosatelliteIndicator + ")";
        }
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
