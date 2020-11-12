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
                        .primaryTumorLocation(Strings.nullToEmpty(patient.baselineData().curatedPrimaryTumor().primaryTumorLocation()))
                        .primaryTumorSubLocation(Strings.nullToEmpty(patient.baselineData()
                                .curatedPrimaryTumor()
                                .primaryTumorSubLocation()))
                        .primaryTumorType(Strings.nullToEmpty(patient.baselineData().curatedPrimaryTumor().primaryTumorType()))
                        .primaryTumorSubType(Strings.nullToEmpty(patient.baselineData().curatedPrimaryTumor().primaryTumorSubType()))
                        .primaryTumorExtraDetails(Strings.nullToEmpty(patient.baselineData()
                                .curatedPrimaryTumor()
                                .primaryTumorExtraDetails()))
                        .doids(extractDoids(patient.baselineData().curatedPrimaryTumor().doidNodes()))
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
}
