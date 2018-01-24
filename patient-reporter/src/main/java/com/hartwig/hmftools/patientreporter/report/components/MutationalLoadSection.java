package com.hartwig.hmftools.patientreporter.report.components;

import java.awt.Color;

import org.jetbrains.annotations.NotNull;

import net.sf.dynamicreports.report.builder.component.ComponentBuilder;

public final class MutationalLoadSection {
    private static final int DRUP_THRESHOLD = 140;

    //MIVO: buffer at edges of slider bar for min/max values. e.g. value of 1 would put min=1 and max=99
    private static final int BUFFER = 3;
    private static final int START = 1;
    private static final int END = 1000;

    @NotNull
    public static ComponentBuilder<?, ?> build(final int mutationalLoad) {
        final int graphValue = computeGraphValue(mutationalLoad);
        final int markerValue = computeGraphValue(DRUP_THRESHOLD);

        final GradientBar gradient =
                ImmutableGradientBar.of(new Color(253, 235, 171), new Color(70, 81, 137), "Low", "High", graphValue, markerValue);
        final SliderSection sliderSection =
                ImmutableSliderSection.of("Tumor Mutational Load", interpret(mutationalLoad), description(), gradient);
        return sliderSection.build();
    }

    @NotNull
    private static String interpret(final int mutationalLoad) {
        if (mutationalLoad > DRUP_THRESHOLD) {
            return "High (" + mutationalLoad + ")";
        } else {
            return "Low (" + mutationalLoad + ")";
        }
    }

    private static int computeGraphValue(final double value) {
        final double logStart = Math.log10(START);
        final double logEnd = Math.log10(END);
        final double logValue = Math.log10(value);
        final double logScaleValue = Math.min(logEnd, Math.max(logStart, logValue));
        final double logIntervalLength = logEnd - logStart;

        final double graphIntervalLength = 100 - 2 * BUFFER;
        return (int) Math.round((logScaleValue - logStart) * graphIntervalLength / logIntervalLength + BUFFER);
    }

    @NotNull
    private static String description() {
        return "Patients with a mutational load over " + DRUP_THRESHOLD + " could be eligible for immunotherapy within the DRUP study. ";
    }
}
