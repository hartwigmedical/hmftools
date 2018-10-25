package com.hartwig.hmftools.patientreporter.report.components;

import java.awt.Color;
import java.text.DecimalFormat;

import org.jetbrains.annotations.NotNull;

import net.sf.dynamicreports.report.builder.component.ComponentBuilder;

public class ChordSection {

    private static final int BUFFER = 3;
    private static final double START = 1E-2;
    private static final double END = 120;

    @NotNull
    public static ComponentBuilder<?, ?> build(final double chordValue) {
        final int graphValue = computeGraphValue(chordValue);

        final GradientBar gradient = ImmutableGradientBar.of(new Color(214, 234, 248), new Color(174, 214, 241), "Low", "High", graphValue);
        final SliderSection sliderSection = ImmutableSliderSection.of("Chord Value",
                interpret(chordValue),
                description(),
                gradient);

        return sliderSection.build();
    }

    @NotNull
    private static String interpret(final double chordValue) {
        return new DecimalFormat("#.####").format(chordValue)  + " value of HR detection";
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
        return Math.log10(value);
    }

    @NotNull
    private static String description() {
        return "";
    }
}
