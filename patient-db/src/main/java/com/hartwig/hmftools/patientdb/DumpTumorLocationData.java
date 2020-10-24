package com.hartwig.hmftools.patientdb;

import java.io.IOException;
import java.util.Collection;
import java.util.List;
import java.util.stream.Collectors;

import com.google.common.base.Strings;
import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.clinical.ImmutablePatientTumorLocation;
import com.hartwig.hmftools.common.clinical.ImmutablePatientTumorLocationV2;
import com.hartwig.hmftools.common.clinical.PatientTumorLocation;
import com.hartwig.hmftools.common.clinical.PatientTumorLocationFile;
import com.hartwig.hmftools.common.clinical.PatientTumorLocationV2;
import com.hartwig.hmftools.common.clinical.PatientTumorLocationV2File;
import com.hartwig.hmftools.common.doid.Doid;
import com.hartwig.hmftools.patientdb.data.Patient;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.jetbrains.annotations.NotNull;

final class DumpTumorLocationData {

    private static final Logger LOGGER = LogManager.getLogger(DumpTumorLocationData.class);

    private DumpTumorLocationData() {
    }

    static void writeCuratedTumorLocationsToTSV(@NotNull String outputFile, @NotNull Collection<Patient> patients) throws IOException {
        List<PatientTumorLocation> tumorLocations = toPatientTumorLocations(patients);
        PatientTumorLocationFile.writeRecordsTSV(outputFile, tumorLocations);
        LOGGER.info(" Written {} tumor locations to {}", tumorLocations.size(), outputFile);
    }

    static void writeCuratedTumorLocationsV2ToTSV(@NotNull String outputFile, @NotNull Collection<Patient> patients) throws IOException {
        List<PatientTumorLocationV2> tumorLocations = toPatientTumorLocationsV2(patients);
        PatientTumorLocationV2File.write(outputFile, tumorLocations);
        LOGGER.info(" Written {} v2 tumor locations to {}", tumorLocations.size(), outputFile);
    }

    @NotNull
    private static List<PatientTumorLocationV2> toPatientTumorLocationsV2(@NotNull Collection<Patient> patients) {
        return patients.stream()
                .map(patient -> ImmutablePatientTumorLocationV2.builder()
                        .patientIdentifier(patient.patientIdentifier())
                        .primaryTumorLocation(Strings.nullToEmpty(patient.baselineData().curatedTumorLocationV2().primaryTumorLocation()))
                        .primaryTumorSubLocation(Strings.nullToEmpty(patient.baselineData()
                                .curatedTumorLocationV2()
                                .primaryTumorSubLocation()))
                        .primaryTumorType(Strings.nullToEmpty(patient.baselineData().curatedTumorLocationV2().primaryTumorType()))
                        .primaryTumorSubType(Strings.nullToEmpty(patient.baselineData().curatedTumorLocationV2().primaryTumorSubType()))
                        .primaryTumorExtraDetails(Strings.nullToEmpty(patient.baselineData()
                                .curatedTumorLocationV2()
                                .primaryTumorExtraDetails()))
                        .doids(extractDoids(patient.baselineData().curatedTumorLocationV2().doids()))
                        .isOverridden(false)
                        .build()).collect(Collectors.toList());
    }

    @NotNull
    private static List<String> extractDoids(@NotNull List<Doid> doidEntries) {
        List<String> doids = Lists.newArrayList();
        for (Doid doidEntry : doidEntries) {
            doids.add(doidEntry.doid());
        }
        return doids;
    }

    @NotNull
    private static List<PatientTumorLocation> toPatientTumorLocations(@NotNull Collection<Patient> patients) {
        return patients.stream()
                .map(patient -> ImmutablePatientTumorLocation.builder()
                        .patientIdentifier(patient.patientIdentifier())
                        .primaryTumorLocation(Strings.nullToEmpty(patient.baselineData().curatedTumorLocation().primaryTumorLocation()))
                        .cancerSubtype(Strings.nullToEmpty(patient.baselineData().curatedTumorLocation().subType()))
                        .build())
                .collect(Collectors.toList());
    }
}
