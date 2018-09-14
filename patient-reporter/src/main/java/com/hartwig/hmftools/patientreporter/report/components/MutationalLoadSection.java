package com.hartwig.hmftools.patientreporter.report.components;

import java.awt.Color;

import com.hartwig.hmftools.common.purple.purity.FittedPurityStatus;
import com.hartwig.hmftools.patientreporter.report.util.PatientReportFormat;

import org.jetbrains.annotations.NotNull;

import net.sf.dynamicreports.report.builder.component.ComponentBuilder;

public final class MutationalLoadSection {
    private static final int DRUP_THRESHOLD = 140;

    //MIVO: buffer at edges of slider bar for min/max values. e.g. value of 1 would put min=1 and max=99
    private static final int BUFFER = 3;
    private static final int START = 1;
    private static final int END = 1000;

    @NotNull
    public static ComponentBuilder<?, ?> build(final int mutationalLoad, @NotNull FittedPurityStatus fitStatus) {
        final int graphValue = computeGraphValue(mutationalLoad);
        final int markerValue = computeGraphValue(DRUP_THRESHOLD);
        final GradientBar gradient =
                ImmutableGradientBar.of(new Color(239, 229, 203), new Color(159, 163, 193), "Low", "High", graphValue, markerValue);
        final SliderSection sliderSection =
                ImmutableSliderSection.of("Tumor Mutational Load", interpret(mutationalLoad, fitStatus), description(), gradient);
        return sliderSection.build();
    }

    @NotNull
    private static String interpret(final int mutationalLoad, @NotNull FittedPurityStatus fitStatus) {
        final String formattedMutationalLoad = PatientReportFormat.correctValueForFitStatus(fitStatus, Integer.toString(mutationalLoad));
        if (mutationalLoad > DRUP_THRESHOLD) {
            return "High (" + formattedMutationalLoad + ")";
        } else {
            return "Low (" + formattedMutationalLoad + ")";
        }
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
        return "The tumor mutational load represents the total number of somatic missense somaticVariants across"
                + " the whole genome of the tumor. Patients with a mutational load over " + DRUP_THRESHOLD
                + " could be eligible for immunotherapy within the DRUP study. ";
    }
}
