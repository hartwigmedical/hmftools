package com.hartwig.hmftools.patientdb;

import java.io.IOException;
import java.util.Collection;
import java.util.List;
import java.util.stream.Collectors;

import com.google.common.base.Strings;
import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.clinical.ImmutablePatientTumorLocation;
import com.hartwig.hmftools.common.clinical.PatientTumorLocation;
import com.hartwig.hmftools.common.clinical.PatientTumorLocationFile;
import com.hartwig.hmftools.common.doid.DoidNode;
import com.hartwig.hmftools.patientdb.data.Patient;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

final class DumpTumorLocationData {

    private static final Logger LOGGER = LogManager.getLogger(DumpTumorLocationData.class);

    private DumpTumorLocationData() {
    }

    static void writeCuratedTumorLocationsToTSV(@NotNull String outputFile, @NotNull Collection<Patient> patients) throws IOException {
        List<PatientTumorLocation> tumorLocations = toPatientTumorLocations(patients);
        PatientTumorLocationFile.write(outputFile, tumorLocations);
        LOGGER.info(" Written {} v2 tumor locations to {}", tumorLocations.size(), outputFile);
    }

    @NotNull
    private static List<PatientTumorLocation> toPatientTumorLocations(@NotNull Collection<Patient> patients) {
        return patients.stream()
                .map(patient -> ImmutablePatientTumorLocation.builder()
                        .patientIdentifier(patient.patientIdentifier())
                        .primaryTumorLocation(Strings.nullToEmpty(patient.baselineData().curatedTumorLocation().primaryTumorLocation()))
                        .primaryTumorSubLocation(Strings.nullToEmpty(patient.baselineData()
                                .curatedTumorLocation()
                                .primaryTumorSubLocation()))
                        .primaryTumorType(Strings.nullToEmpty(patient.baselineData().curatedTumorLocation().primaryTumorType()))
                        .primaryTumorSubType(Strings.nullToEmpty(patient.baselineData().curatedTumorLocation().primaryTumorSubType()))
                        .primaryTumorExtraDetails(Strings.nullToEmpty(patient.baselineData()
                                .curatedTumorLocation()
                                .primaryTumorExtraDetails()))
                        .doids(extractDoids(patient.baselineData().curatedTumorLocation().doidNodes()))
                        .isOverridden(false)
                        .build()).collect(Collectors.toList());
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
