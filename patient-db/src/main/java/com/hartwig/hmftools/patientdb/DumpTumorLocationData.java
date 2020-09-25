package com.hartwig.hmftools.patientdb;

import java.io.File;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
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

    static void writeCuratedTumorLocationsToCSV(@NotNull String tumorLocationOutputDir,
            @NotNull Collection<Patient> patients) throws IOException {
        String outputFile = fileLocation(tumorLocationOutputDir, "curatedTumorLocations.csv");
        List<PatientTumorLocation> tumorLocations = patients.stream()
                .map(patient -> ImmutablePatientTumorLocation.of(patient.patientIdentifier(),
                        Strings.nullToEmpty(patient.baselineData().curatedTumorLocation().primaryTumorLocation()),
                        Strings.nullToEmpty(patient.baselineData().curatedTumorLocation().subType())))
                .collect(Collectors.toList());

        PatientTumorLocation.writeRecords(outputFile, tumorLocations);
        LOGGER.info(" Written {} tumor locations to {}.", tumorLocations.size(), outputFile);
    }

    @NotNull
    private static String fileLocation(@NotNull String outputDir, @NotNull String fileName) {
        return outputDir + File.separator + fileName;
    }

    private static void updateSymlink(@NotNull String linkName, @NotNull String fileName) {
        Path linkPath = Paths.get(linkName);
        try {
            Files.deleteIfExists(linkPath);
            Files.createSymbolicLink(linkPath, Paths.get(fileName));
        } catch (IOException e) {
            LOGGER.warn("Failed to update symlink {}. Cause: {}.", linkName, e.getMessage());
        }
    }
}
