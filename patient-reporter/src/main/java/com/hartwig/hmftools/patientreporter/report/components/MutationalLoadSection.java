package com.hartwig.hmftools.patientreporter.report.components;

import java.awt.Color;

import org.jetbrains.annotations.NotNull;

import net.sf.dynamicreports.report.builder.component.ComponentBuilder;

public final class MutationalLoadSection {
    private static final double THRESHOLD = 140;
    //MIVO: buffer at edges of slider bar for min/max values. e.g. value of 1 would put min=1 and max=99
    private static final int BUFFER = 3;
    private static final int START = 0;
    private static final int END = 3;

    @NotNull
    public static ComponentBuilder<?, ?> build(final int mutationalLoad) {
        final int graphValue = computeGraphValue(mutationalLoad);
        final int markerValue = computeGraphValue(THRESHOLD);
        final GradientBar gradient =
                ImmutableGradientBar.of(new Color(253, 235, 171), new Color(70, 81, 137), "Normal", "High", graphValue, markerValue);
        final SliderSection sliderSection =
                ImmutableSliderSection.of("Mutational Load", interpret(mutationalLoad), description(), gradient);
        return sliderSection.build();
    }

    @NotNull
    private static String interpret(final int mutationalLoad) {
        if (mutationalLoad > THRESHOLD) {
            return "High (" + mutationalLoad + ")";
        } else {
            return "Normal (" + mutationalLoad + ")";
        }
    }

    private static int computeGraphValue(final double value) {
        final double logValue = Math.log10(value);
        final double logScaleValue = Math.min(END, Math.max(START, logValue));
        final int logIntervalLength = END - START;
        final int graphIntervalLength = 100 - 2 * BUFFER;
        return (int) Math.round((logScaleValue - START) * graphIntervalLength / logIntervalLength + BUFFER);
    }

    @NotNull
    private static String description() {
        return "Mutational load represents the total number of observed somatic variants across the entire protein coding region of the genome. "
                + "Somatic variants only include single nucleotide variants (excluding insertion/deletions) and represent potential neoantigens. "
                + "Patients with a mutational load over " + THRESHOLD + " could be eligible for immunotherapy within the DRUP study. ";
        //                + "The arrow on the bar scale indicates the relative position of the analyzed sample compared to all samples (upper bar) "
        //                + "and compared to all ‘on-tumor-type’ samples (lower bar) in the database.";
    }
}
