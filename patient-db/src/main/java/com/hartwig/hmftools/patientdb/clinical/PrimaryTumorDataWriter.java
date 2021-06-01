package com.hartwig.hmftools.patientdb.clinical;

import java.io.IOException;
import java.util.Collection;
import java.util.List;
import java.util.Map;
import java.util.stream.Collectors;

import com.google.common.base.Strings;
import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import com.hartwig.hmftools.common.clinical.ImmutablePatientPrimaryTumor;
import com.hartwig.hmftools.common.clinical.PatientPrimaryTumor;
import com.hartwig.hmftools.common.clinical.PatientPrimaryTumorFile;
import com.hartwig.hmftools.common.doid.DoidNode;
import com.hartwig.hmftools.common.lims.reportingdb.ReportingDatabase;
import com.hartwig.hmftools.common.lims.reportingdb.ReportingEntry;
import com.hartwig.hmftools.patientdb.clinical.curators.PatientTumorCurationStatus;
import com.hartwig.hmftools.patientdb.clinical.curators.PatientTumorCurationStatusFile;
import com.hartwig.hmftools.patientdb.clinical.datamodel.Patient;
import com.hartwig.hmftools.patientdb.clinical.datamodel.SampleData;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

public final class PrimaryTumorDataWriter {

    private static final Logger LOGGER = LogManager.getLogger(PrimaryTumorDataWriter.class);

    private PrimaryTumorDataWriter() {
    }

    public static void write(@NotNull ClinicalAlgoConfig config, @NotNull Map<String, List<SampleData>> samplesPerPatient,
            @NotNull List<Patient> patients) throws IOException {
        LOGGER.info("Check for missing curation tumor location when info is known");
        Map<String, PatientTumorCurationStatus> patientTumorCurationStatusMap =
                generatePatientTumorCurationStatusMap(samplesPerPatient, patients, config.reportingDbTsv());

        PrimaryTumorDataWriter.writePatientTumorCurationStatesToTSV(config.patientTumorCurationStatusTsv(), patientTumorCurationStatusMap);
        PrimaryTumorDataWriter.writeCuratedPrimaryTumorsToTSV(config.curatedPrimaryTumorTsv(), patients);
    }

    @NotNull
    private static Map<String, PatientTumorCurationStatus> generatePatientTumorCurationStatusMap(
            @NotNull Map<String, List<SampleData>> samplesPerPatient, @NotNull List<Patient> patients, @NotNull String reportingDbTsv)
            throws IOException {
        Map<String, PatientTumorCurationStatus> patientTumorCurationStatusMap = Maps.newHashMap();

        List<String> reportedBarcodes = Lists.newArrayList();
        for (ReportingEntry entry : ReportingDatabase.read(reportingDbTsv)) {
            reportedBarcodes.add(entry.tumorBarcode());
        }

        List<String> patientsWithSamplesToBeReported = Lists.newArrayList();
        for (Map.Entry<String, List<SampleData>> entry : samplesPerPatient.entrySet()) {
            String patientId = entry.getKey();
            for (SampleData sample : entry.getValue()) {
                if (!sample.requiresCuratedPrimaryTumor()) {
                    patientTumorCurationStatusMap.put(patientId, PatientTumorCurationStatus.NEEDS_NO_CURATED_PRIMARY_TUMOR);
                } else if (sample.isSomaticTumorSample()) {
                    if (reportedBarcodes.contains(sample.sampleBarcode())) {
                        if (allowExceptionsForPatient(patientId)) {
                            patientTumorCurationStatusMap.put(patientId, PatientTumorCurationStatus.ALREADY_REPORTED);
                            patientsWithSamplesToBeReported.add(patientId);
                        } else {
                            patientsWithSamplesToBeReported.add(patientId);
                        }
                    } else {
                        patientsWithSamplesToBeReported.add(patientId);
                    }
                }
            }
        }

        for (String patientId : patientsWithSamplesToBeReported) {
            Patient patient = findByPatientId(patients, patientId);
            if (patient != null) {
                String tumorLocationSearchTerm = patient.baselineData().curatedPrimaryTumor().searchTerm();
                if (tumorLocationSearchTerm != null && !tumorLocationSearchTerm.isEmpty()) {
                    if (patient.baselineData().curatedPrimaryTumor().location() == null) {
                        LOGGER.warn("Could not curate patient {} for primary tumor '{}'",
                                patient.patientIdentifier(),
                                tumorLocationSearchTerm);
                    }
                } else {
                    if (!allowExceptionsForPatient(patient.patientIdentifier())) {
                        LOGGER.warn("Could not find input tumor location for patient {}", patient.patientIdentifier());
                    } else {
                        patientTumorCurationStatusMap.put(patient.patientIdentifier(), PatientTumorCurationStatus.MISSING_TUMOR_CURATION);
                    }
                }
            } else {
                if (allowExceptionsForPatient(patientId)) {
                    patientTumorCurationStatusMap.put(patientId, PatientTumorCurationStatus.NOT_RESOLVED);
                }
            }
        }

        return patientTumorCurationStatusMap;
    }

    @Nullable
    private static Patient findByPatientId(@NotNull List<Patient> patients, @NotNull String patientId) {
        for (Patient patient : patients) {
            if (patient.patientIdentifier().equals(patientId)) {
                return patient;
            }
        }

        return null;
    }

    private static boolean allowExceptionsForPatient(@NotNull String patientId) {
        return !patientId.startsWith("CORE") && !patientId.startsWith("WIDE");
    }

    private static void writePatientTumorCurationStatesToTSV(@NotNull String outputTsv,
            @NotNull Map<String, PatientTumorCurationStatus> patientTumorCurationStatusMap) throws IOException {
        LOGGER.info("Writing patient tumor curation status");
        PatientTumorCurationStatusFile.write(outputTsv, patientTumorCurationStatusMap);
        LOGGER.info(" Written {} patient tumor curation states to {}", patientTumorCurationStatusMap.size(), outputTsv);
    }

    private static void writeCuratedPrimaryTumorsToTSV(@NotNull String outputTsv, @NotNull Collection<Patient> patients)
            throws IOException {
        List<PatientPrimaryTumor> primaryTumors = toPatientPrimaryTumors(patients);
        LOGGER.info("Writing curated primary tumors");
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
