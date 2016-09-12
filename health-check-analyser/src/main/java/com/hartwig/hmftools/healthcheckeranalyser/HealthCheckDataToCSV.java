package com.hartwig.hmftools.healthcheckeranalyser;

import java.util.Set;
import java.util.SortedSet;

import com.google.common.collect.Sets;
import com.hartwig.hmftools.healthcheckeranalyser.model.HealthCheckReport;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.jetbrains.annotations.NotNull;

final class HealthCheckDataToCSV {

    private static final Logger LOGGER = LogManager.getLogger(HealthCheckDataToCSV.class);

    private HealthCheckDataToCSV() {
    }

    @NotNull
    static String header(@NotNull HealthCheckReport report) {
        // KODU: We assume every report generates the same header.
        String header = "SAMPLE,RUN_DATE,PIPELINE_VERSION";
        for (String check : getSampleChecks(report)) {
            header += ("," + check);
        }

        for (String check : getPatientChecks(report)) {
            header += ("," + check);
        }

        return header;
    }

    @NotNull
    static String refSample(@NotNull HealthCheckReport report) {
        String refSample = report.refSample() + "," + report.runDate() + "," + report.pipelineVersion();
        for (String check : getSampleChecks(report)) {
            refSample += ("," + report.refChecks().get(check));
        }
        for (String ignored : getPatientChecks(report)) {
            refSample += (",-");
        }

        return refSample;
    }

    @NotNull
    static String tumorSample(@NotNull HealthCheckReport report) {
        String tumorSample = report.tumorSample() + "," + report.runDate() + "," + report.pipelineVersion();
        for (String check : getSampleChecks(report)) {
            tumorSample += ("," + report.tumorChecks().get(check));
        }
        for (String check : getPatientChecks(report)) {
            tumorSample += ("," + report.patientChecks().get(check));
        }

        return tumorSample;
    }

    @NotNull
    private static Set<String> getSampleChecks(@NotNull HealthCheckReport report) {
        SortedSet<String> refChecks = Sets.newTreeSet(report.refChecks().keySet());
        SortedSet<String> tumorChecks = Sets.newTreeSet(report.tumorChecks().keySet());
        if (!refChecks.equals(tumorChecks)) {
            LOGGER.warn("RefChecks and TumorChecks do not match!");
        }
        return refChecks;
    }

    @NotNull
    private static Set<String> getPatientChecks(@NotNull HealthCheckReport report) {
        return Sets.newTreeSet(report.patientChecks().keySet());
    }
}
