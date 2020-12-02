package com.hartwig.hmftools.patientdb;

import java.io.IOException;
import java.util.Collection;
import java.util.List;
import java.util.stream.Collectors;

import com.google.common.base.Strings;
import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.clinical.ImmutablePatientPrimaryTumor;
import com.hartwig.hmftools.common.clinical.PatientPrimaryTumor;
import com.hartwig.hmftools.common.clinical.PatientPrimaryTumorFile;
import com.hartwig.hmftools.common.doid.DoidNode;
import com.hartwig.hmftools.patientdb.data.Patient;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

final class DumpPrimaryTumorData {

    private static final Logger LOGGER = LogManager.getLogger(DumpPrimaryTumorData.class);

    private DumpPrimaryTumorData() {
    }

    static void writeCuratedPrimaryTumorsToTSV(@NotNull String outputTsv, @NotNull Collection<Patient> patients) throws IOException {
        List<PatientPrimaryTumor> primaryTumors = toPatientPrimaryTumors(patients);
        PatientPrimaryTumorFile.write(outputTsv, primaryTumors);
        LOGGER.info(" Written {} primary tumors to {}", primaryTumors.size(), outputTsv);
    }

    @NotNull
    private static List<PatientPrimaryTumor> toPatientPrimaryTumors(@NotNull Collection<Patient> patients) {
        return patients.stream()
                .map(patient -> ImmutablePatientPrimaryTumor.builder()
                        .patientIdentifier(patient.patientIdentifier())
                        .location(Strings.nullToEmpty(patient.baselineData().curatedPrimaryTumor().location()))
                        .subLocation(Strings.nullToEmpty(patient.baselineData().curatedPrimaryTumor().subLocation()))
                        .type(Strings.nullToEmpty(patient.baselineData().curatedPrimaryTumor().type()))
                        .subType(Strings.nullToEmpty(patient.baselineData().curatedPrimaryTumor().subType()))
                        .extraDetails(Strings.nullToEmpty(patient.baselineData().curatedPrimaryTumor().extraDetails()))
                        .doids(extractDoids(patient.baselineData().curatedPrimaryTumor().doidNodes()))
                        .snomedConceptIds(nullToEmpty(patient.baselineData().curatedPrimaryTumor().snomedConceptIds()))
                        .isOverridden(false)
                        .build())
                .collect(Collectors.toList());
    }

    @NotNull
    private static List<String> extractDoids(@Nullable List<DoidNode> doidNodes) {
        if (doidNodes == null) {
            return Lists.newArrayList();
        }

        List<String> doids = Lists.newArrayList();
        for (DoidNode doidNode : doidNodes) {
            doids.add(doidNode.doid());
        }
        return doids;
    }

    @NotNull
    private static List<String> nullToEmpty(@Nullable List<String> strings) {
        return strings != null ? strings : Lists.newArrayList();
    }
}
