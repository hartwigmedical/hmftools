package com.hartwig.hmftools.healthcheckreader.model;

import java.io.FileNotFoundException;
import java.io.FileReader;
import java.util.Map;

import com.google.common.collect.Maps;
import com.google.gson.Gson;
import com.google.gson.GsonBuilder;
import com.google.gson.JsonArray;
import com.google.gson.JsonElement;
import com.google.gson.JsonObject;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.apache.logging.log4j.util.Strings;
import org.jetbrains.annotations.NotNull;

public final class HealthCheckReportFactory {

    private static final Logger LOGGER = LogManager.getLogger(HealthCheckReportFactory.class);
    private static final Gson GSON = new GsonBuilder().create();

    // KODU: The format of health checks depends on github health-checks project.
    private static final String HEALTH_CHECK_MAIN_OBJECT = "health_checks";

    private static final String PATIENT_CHECKS_IDENTIFIER_1 = "check";
    private static final String PATIENT_CHECKS_IDENTIFIER_2 = "checks";
    private static final String REF_SAMPLE_CHECKS_IDENTIFIER = "ref_sample_checks";
    private static final String TUMOR_SAMPLE_CHECKS_IDENTIFIER = "tumor_sample_checks";

    private static final String CHECK_SAMPLE_ID = "sample_id";
    private static final String CHECK_NAME = "check_name";
    private static final String CHECK_VALUE = "value";

    private HealthCheckReportFactory() {
    }

    @NotNull
    public static HealthCheckReport fromHealthCheckReport(@NotNull final String path) throws FileNotFoundException {
        final JsonObject json = GSON.fromJson(new FileReader(path), JsonObject.class);
        final JsonArray checks = json.get(HEALTH_CHECK_MAIN_OBJECT).getAsJsonArray();

        String refSample = Strings.EMPTY;
        String tumorSample = Strings.EMPTY;
        final Map<String, String> refChecks = Maps.newHashMap();
        final Map<String, String> tumorChecks = Maps.newHashMap();
        final Map<String, String> patientChecks = Maps.newHashMap();

        for (int i = 0; i < checks.size(); i++) {
            JsonObject category = checks.get(i).getAsJsonObject();
            for (final Map.Entry<String, JsonElement> entry : category.entrySet()) {
                final JsonObject values = entry.getValue().getAsJsonObject();
                if (values.has(PATIENT_CHECKS_IDENTIFIER_1)) {
                    patientChecks.putAll(extractChecks(values.get(PATIENT_CHECKS_IDENTIFIER_1)));
                } else if (values.has(PATIENT_CHECKS_IDENTIFIER_2)) {
                    patientChecks.putAll(extractChecks(values.get(PATIENT_CHECKS_IDENTIFIER_2)));
                } else if (values.has(REF_SAMPLE_CHECKS_IDENTIFIER) && values.has(TUMOR_SAMPLE_CHECKS_IDENTIFIER)) {
                    if (refSample.equals(Strings.EMPTY)) {
                        refSample = extractSampleId(values.get(REF_SAMPLE_CHECKS_IDENTIFIER));
                    }
                    refChecks.putAll(extractChecks(values.get(REF_SAMPLE_CHECKS_IDENTIFIER)));
                    if (tumorSample.equals(Strings.EMPTY)) {
                        tumorSample = extractSampleId(values.get(TUMOR_SAMPLE_CHECKS_IDENTIFIER));
                    }
                    tumorChecks.putAll(extractChecks(values.get(TUMOR_SAMPLE_CHECKS_IDENTIFIER)));
                } else {
                    LOGGER.error("Unrecognized category: " + values);
                }
            }
        }

        if (refChecks.size() != tumorChecks.size()) {
            LOGGER.error("Size of ref checks and tumor checks do not match");
        }

        return new HealthCheckReport(refSample, tumorSample, refChecks, tumorChecks, patientChecks);
    }

    @NotNull
    private static Map<String, String> extractChecks(@NotNull final JsonElement element) {
        final Map<String, String> elements = Maps.newHashMap();
        if (element.isJsonObject()) {
            final JsonObject object = (JsonObject) element;
            elements.put(object.get(CHECK_NAME).getAsString(), object.get(CHECK_VALUE).getAsString());
        } else if (element.isJsonArray()) {
            final JsonArray array = (JsonArray) element;
            for (int i = 0; i < array.size(); i++) {
                final JsonObject object = (JsonObject) array.get(i);
                elements.put(object.get(CHECK_NAME).getAsString(), object.get(CHECK_VALUE).getAsString());
            }
        }
        return elements;
    }

    @NotNull
    private static String extractSampleId(@NotNull final JsonElement element) {
        final JsonArray array = (JsonArray) element;
        final JsonObject firstObject = (JsonObject) array.get(0);
        return firstObject.get(CHECK_SAMPLE_ID).getAsString();
    }
}
