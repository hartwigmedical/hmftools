package com.hartwig.hmftools.patientdb;

import java.io.File;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.time.LocalDate;
import java.time.format.DateTimeFormatter;
import java.util.Collection;
import java.util.List;
import java.util.Optional;
import java.util.stream.Collectors;

import com.google.common.base.Strings;
import com.hartwig.hmftools.common.ecrf.projections.ImmutablePatientTumorLocation;
import com.hartwig.hmftools.common.ecrf.projections.ImmutablePortalClinicalData;
import com.hartwig.hmftools.common.ecrf.projections.PatientTumorLocation;
import com.hartwig.hmftools.common.ecrf.projections.PortalClinicalData;
import com.hartwig.hmftools.patientdb.data.Patient;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.jetbrains.annotations.NotNull;

final class DumpClinicalData {

    private static final Logger LOGGER = LogManager.getLogger(DumpClinicalData.class);

    private DumpClinicalData() {
    }

    static void writeClinicalDumps(@NotNull final String csvOutputDir, @NotNull final Collection<Patient> patients,
            @NotNull final Optional<String> tumorLocationLink, @NotNull final Optional<String> portalDataLink) throws IOException {
        writeCuratedTumorLocationsToCSV(csvOutputDir, tumorLocationLink, patients);
        writePortalClinicalData(csvOutputDir, portalDataLink, patients);
    }

    private static void writeCuratedTumorLocationsToCSV(@NotNull final String csvOutputDir, @NotNull final Optional<String> linkName,
            @NotNull final Collection<Patient> patients) throws IOException {
        final String outputFile = fileLocation(csvOutputDir, "_curatedTumorLocations.csv");
        LOGGER.info("Writing curated tumor locations to csv in {}.", csvOutputDir);
        final List<PatientTumorLocation> tumorLocations = patients.stream()
                .map(patient -> ImmutablePatientTumorLocation.of(patient.patientIdentifier(),
                        Strings.nullToEmpty(patient.baselineData().curatedTumorLocation().primaryTumorLocation()),
                        Strings.nullToEmpty(patient.baselineData().curatedTumorLocation().subType())))
                .collect(Collectors.toList());
        PatientTumorLocation.writeRecords(outputFile, tumorLocations);
        linkName.ifPresent(link -> updateSymlink(csvOutputDir + File.separator + link, outputFile));
        LOGGER.info("Written {} records to {}.", tumorLocations.size(), outputFile);
    }

    private static void writePortalClinicalData(@NotNull final String csvOutputDir, @NotNull final Optional<String> linkName,
            @NotNull final Collection<Patient> patients) throws IOException {
        final String outputFile = fileLocation(csvOutputDir, "_portal.csv");
        LOGGER.info("Writing portal clinical data to csv in {}.", csvOutputDir);
        final List<PortalClinicalData> portalData = patients.stream()
                .flatMap(patient -> patient.sequencedBiopsies()
                        .stream()
                        .map(sampleData -> ImmutablePortalClinicalData.of(patient.patientIdentifier(),
                                sampleData.sampleId(),
                                patient.baselineData().gender(),
                                patient.baselineData().birthYear(),
                                patient.baselineData().registrationDate(),
                                patient.baselineData().curatedTumorLocation().primaryTumorLocation(),
                                getBiopsyType(patient, sampleData.sampleId()))))
                .collect(Collectors.toList());
        PortalClinicalData.writeRecords(outputFile, portalData);
        linkName.ifPresent(link -> updateSymlink(csvOutputDir + File.separator + link, outputFile));
        LOGGER.info("Written {} sample records to {}.", portalData.size(), outputFile);
    }

    @NotNull
    private static String getBiopsyType(@NotNull final Patient patient, @NotNull final String sampleId) {
        return patient.clinicalBiopsies().stream().filter(biopsyData -> {
            final String clinicalSampleId = biopsyData.sampleId();
            return clinicalSampleId != null && clinicalSampleId.equals(sampleId);
        }).map(biopsyData -> Strings.nullToEmpty(biopsyData.curatedType())).findFirst().orElse("");
    }

    @NotNull
    private static String fileLocation(@NotNull final String outputDir, @NotNull final String suffix) {
        final String fileName = LocalDate.now().format(DateTimeFormatter.ISO_LOCAL_DATE) + suffix;
        return outputDir + File.separator + fileName;
    }

    private static void updateSymlink(@NotNull final String linkName, @NotNull final String fileName) {
        final Path linkPath = Paths.get(linkName);
        try {
            Files.deleteIfExists(linkPath);
            Files.createSymbolicLink(linkPath, Paths.get(fileName));
        } catch (IOException e) {
            LOGGER.warn("Failed to update symlink {}. Cause: {}.", linkName, e.getMessage());
        }
    }
}
