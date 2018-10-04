package com.hartwig.hmftools.patientreporter.report.components;

import java.awt.Color;
import java.text.DecimalFormat;

import com.hartwig.hmftools.common.purple.purity.FittedPurityStatus;
import com.hartwig.hmftools.patientreporter.report.util.PatientReportFormat;

import org.jetbrains.annotations.NotNull;

import net.sf.dynamicreports.report.builder.component.ComponentBuilder;

public class TumorMutationBurdenSection {
    private static final int BUFFER = 3;
    private static final double START = 1E-2;
    private static final double END = 100;
    private static final double TMB_THRESHOLD = 10D;

    @NotNull
    public static ComponentBuilder<?, ?> build(final double tumorMutationalBurdenIndicator, @NotNull FittedPurityStatus fitStatus) {
        final int graphValue = computeGraphValue(tumorMutationalBurdenIndicator);
        final int markerValue = computeGraphValue(TMB_THRESHOLD);

        final GradientBar gradient =
                ImmutableGradientBar.of(new Color(239, 239, 239), new Color(171, 191, 171), "Low", "High", graphValue, markerValue);
        final SliderSection sliderSection = ImmutableSliderSection.of("Tumor Mutational Burden",
                interpret(tumorMutationalBurdenIndicator, fitStatus),
                description(),
                gradient);
        return sliderSection.build();
    }

    @NotNull
    private static String interpret(final double tumorMutationalBurden, @NotNull FittedPurityStatus fitStatus) {
        final String formattedTumorMutationBurden =
                PatientReportFormat.correctValueForFitStatus(fitStatus, new DecimalFormat("#.####").format(tumorMutationalBurden));
        if (tumorMutationalBurden > TMB_THRESHOLD) {
            return "High (" + formattedTumorMutationBurden + ")";
        } else {
            return "Low (" + formattedTumorMutationBurden + ")";
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
        return "The tumor mutational burden score represents the number of all somatic variants "
                + "across the whole genome of the tumor per Mb. Tumors with a score greater than " + TMB_THRESHOLD
                + " are considered as high tumor mutation burden (TMB).";
    }
}
