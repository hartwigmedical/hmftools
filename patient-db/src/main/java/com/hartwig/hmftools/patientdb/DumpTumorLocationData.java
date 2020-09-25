package com.hartwig.hmftools.patientdb;

import java.io.IOException;
import java.util.Collection;
import java.util.List;
import java.util.stream.Collectors;

import com.google.common.base.Strings;
import com.hartwig.hmftools.common.ecrf.projections.ImmutablePatientTumorLocation;
import com.hartwig.hmftools.common.ecrf.projections.PatientTumorLocation;
import com.hartwig.hmftools.patientdb.data.Patient;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.jetbrains.annotations.NotNull;

final class DumpTumorLocationData {

    private static final Logger LOGGER = LogManager.getLogger(DumpTumorLocationData.class);

    private DumpTumorLocationData() {
    }

    static void writeCuratedTumorLocationsToCSV(@NotNull String outputFile, @NotNull Collection<Patient> patients) throws IOException {
        List<PatientTumorLocation> tumorLocations = patients.stream()
                .map(patient -> ImmutablePatientTumorLocation.of(patient.patientIdentifier(),
                        Strings.nullToEmpty(patient.baselineData().curatedTumorLocation().primaryTumorLocation()),
                        Strings.nullToEmpty(patient.baselineData().curatedTumorLocation().subType())))
                .collect(Collectors.toList());

        PatientTumorLocation.writeRecordsCSV(outputFile, tumorLocations);
        LOGGER.info(" Written {} tumor locations to {}.", tumorLocations.size(), outputFile);
    }

    static void writeCuratedTumorLocationsToTSV(@NotNull String outputFile, @NotNull Collection<Patient> patients) throws IOException {
        List<PatientTumorLocation> tumorLocations = patients.stream()
                .map(patient -> ImmutablePatientTumorLocation.of(patient.patientIdentifier(),
                        Strings.nullToEmpty(patient.baselineData().curatedTumorLocation().primaryTumorLocation()),
                        Strings.nullToEmpty(patient.baselineData().curatedTumorLocation().subType())))
                .collect(Collectors.toList());

        PatientTumorLocation.writeRecordsTSV(outputFile, tumorLocations);
        LOGGER.info(" Written {} tumor locations to {}.", tumorLocations.size(), outputFile);
    }
}
