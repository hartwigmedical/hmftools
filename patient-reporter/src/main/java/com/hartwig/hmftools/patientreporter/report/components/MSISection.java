package com.hartwig.hmftools.patientreporter.report.components;

import java.awt.Color;
import java.text.DecimalFormat;

import org.jetbrains.annotations.NotNull;

import net.sf.dynamicreports.report.builder.component.ComponentBuilder;

public class MSISection {
    private static final double THRESHOLD = 0.909;
    //MIVO: buffer at edges of slider bar for min/max values. e.g. value of 1 would put min=1 and max=99
    private static final int BUFFER = 3;
    private static final int START = -4;
    private static final int END = 3;

    @NotNull
    public static ComponentBuilder<?, ?> build(final double msiValue) {
        final int graphValue = computeGraphValue(msiValue);
        final int markerValue = computeGraphValue(THRESHOLD);
        final GradientBar gradient = ImmutableGradientBar.of(Color.YELLOW, Color.RED, "MSS", "MSI-H", graphValue, markerValue);
        final SliderSection sliderSection =
                ImmutableSliderSection.of("Microsatellite Instability", interpretMSI(msiValue), description(), gradient);
        return sliderSection.build();
    }

    @NotNull
    private static String interpretMSI(final double msiValue) {
        final String formattedMsi = new DecimalFormat("#.####").format(msiValue);
        if (msiValue > THRESHOLD) {
            return "Unstable (" + formattedMsi + ")";
        } else {
            return "Stable (" + formattedMsi + ")";
        }
    }

    private static int computeGraphValue(final double value) {
        final double logValue = Math.log10(value);
        final double logScaleValue = Math.min(END, Math.max(START, logValue));
        final int logIntervalLength = END - START;
        final int graphIntervalLength = 100 - 2 * BUFFER;
        return (int) Math.round((logScaleValue - START) * graphIntervalLength / logIntervalLength + BUFFER);
    }

    private static String description() {
        return "Microsatellite instability scores are calculated using WGS data and based on the MSIseq classifier as described in "
                + "Huang MN et al, Sci Rep. 5: 13321 (2015). The arrow on the bar scale indicates the relative position of the "
                + "analyzed sample compared to all samples in the  database. Samples with a score greater than " + THRESHOLD
                + " are considered microsatellite instable (MSI-H).";
    }

}
