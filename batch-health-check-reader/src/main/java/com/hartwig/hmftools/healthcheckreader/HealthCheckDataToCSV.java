package com.hartwig.hmftools.healthcheckreader;

import java.util.Set;
import java.util.SortedSet;

import com.google.common.collect.Sets;
import com.hartwig.hmftools.healthcheckreader.model.HealthCheckReport;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.jetbrains.annotations.NotNull;

final class HealthCheckDataToCSV {

    private static final Logger LOGGER = LogManager.getLogger(HealthCheckDataToCSV.class);

    private HealthCheckDataToCSV() {
    }

    @NotNull
    static String header(@NotNull final HealthCheckReport report) {
        // KODU: We assume every report generates the same header.
        String header = "SAMPLE";
        for (final String check : getSampleChecks(report)) {
            header += ("," + check);
        }

        for (final String check : getPatientChecks(report)) {
            header += ("," + check);
        }

        return header;
    }

    @NotNull
    static String refSample(@NotNull final HealthCheckReport report) {
        String refSample = report.refSample();
        for (final String check : getSampleChecks(report)) {
            refSample += ("," + report.refChecks().get(check));
        }
        for (final String ignored : getPatientChecks(report)) {
            refSample += (",-");
        }

        return refSample;
    }

    @NotNull
    static String tumorSample(@NotNull final HealthCheckReport report) {
        String tumorSample = report.tumorSample();
        for (final String check : getSampleChecks(report)) {
            tumorSample += ("," + report.tumorChecks().get(check));
        }
        for (final String check : getPatientChecks(report)) {
            tumorSample += ("," + report.patientChecks().get(check));
        }

        return tumorSample;
    }

    @NotNull
    private static Set<String> getSampleChecks(@NotNull final HealthCheckReport report) {
        final Set<String> refChecks = Sets.newTreeSet(report.refChecks().keySet());
        final Set<String> tumorChecks = Sets.newTreeSet(report.tumorChecks().keySet());
        if (!refChecks.equals(tumorChecks)) {
            LOGGER.warn("RefChecks and TumorChecks do not match!");
        }
        return refChecks;
    }

    @NotNull
    private static Set<String> getPatientChecks(@NotNull final HealthCheckReport report) {
        return Sets.newTreeSet(report.patientChecks().keySet());
    }
}
