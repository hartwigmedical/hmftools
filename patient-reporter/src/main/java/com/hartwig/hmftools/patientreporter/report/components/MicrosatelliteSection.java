package com.hartwig.hmftools.patientreporter.report.components;

import java.awt.Color;
import java.text.DecimalFormat;

import org.jetbrains.annotations.NotNull;

import net.sf.dynamicreports.report.builder.component.ComponentBuilder;

public final class MicrosatelliteSection {
    private static final double MSI_THRESHOLD = 4D;

    //MIVO: buffer at edges of slider bar for min/max values. e.g. value of 1 would put min=1 and max=99
    private static final int BUFFER = 3;
    private static final double START = 1E-2;
    private static final double END = 100;

    @NotNull
    public static ComponentBuilder<?, ?> build(final double microsatelliteIndicator) {
        final int graphValue = computeGraphValue(microsatelliteIndicator);
        final int markerValue = computeGraphValue(MSI_THRESHOLD);

        final GradientBar gradient = ImmutableGradientBar.of(Color.YELLOW, Color.RED, "MSS", "MSI", graphValue, markerValue);
        final SliderSection sliderSection = ImmutableSliderSection.of("Microsatellite Status",
                interpretMicrosatelliteStatus(microsatelliteIndicator),
                description(),
                gradient);
        return sliderSection.build();
    }

    @NotNull
    private static String interpretMicrosatelliteStatus(final double microsatelliteIndicator) {
        final String formattedMsi = new DecimalFormat("#.####").format(microsatelliteIndicator);
        if (microsatelliteIndicator > MSI_THRESHOLD) {
            return "Unstable (" + formattedMsi + ")";
        } else {
            return "Stable (" + formattedMsi + ")";
        }
    }

    private static int computeGraphValue(final double value) {
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
        return "Microsatellite instability scores are calculated using WGS data and based on the MSIseq classifier as described in "
                + "Huang MN et al, Sci Rep. 5: 13321 (2015). The arrow on the bar scale indicates the relative position of the "
                + "analyzed sample compared to all samples in the  database. Samples with a score greater than " + MSI_THRESHOLD
                + " are considered microsatellite unstable (MSI-H).";
    }
}
