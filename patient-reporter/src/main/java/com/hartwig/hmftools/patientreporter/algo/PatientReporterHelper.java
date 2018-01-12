package com.hartwig.hmftools.patientreporter.algo;

import java.io.File;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.List;
import java.util.Optional;

import com.hartwig.hmftools.common.ecrf.CpctEcrfModel;
import com.hartwig.hmftools.common.ecrf.datamodel.EcrfPatient;
import com.hartwig.hmftools.common.exception.HartwigException;
import com.hartwig.hmftools.common.gene.GeneCopyNumber;
import com.hartwig.hmftools.common.gene.GeneCopyNumberFile;
import com.hartwig.hmftools.common.purple.copynumber.PurpleCopyNumber;
import com.hartwig.hmftools.common.purple.copynumber.PurpleCopyNumberFile;
import com.hartwig.hmftools.common.purple.purity.FittedPurityFile;
import com.hartwig.hmftools.common.purple.purity.PurityContext;
import com.hartwig.hmftools.common.variant.SomaticVariant;
import com.hartwig.hmftools.common.variant.SomaticVariantFactory;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.apache.logging.log4j.util.Strings;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

final class PatientReporterHelper {

    private static final Logger LOGGER = LogManager.getLogger(PatientReporterHelper.class);

    private static final String SOMATIC_SNV_EXTENSION = "_post_processed_v2.1.vcf.gz";
    private static final String PURPLE_DIRECTORY = "purple";
    private static final String SV_EXTENSION = "_somaticSV_bpi.vcf";

    private static final String TUMOR_TYPE_ECRF_FIELD = "BASELINE.CARCINOMA.CARCINOMA.PTUMLOC";
    private static final String TUMOR_TYPE_OTHER_ECRF_FIELD = "BASELINE.CARCINOMA.CARCINOMA.PTUMLOCS";

    private PatientReporterHelper() {
    }

    @NotNull
    static PurityContext loadPurity(@NotNull final String runDirectory, @NotNull final String sample) throws IOException {
        final String cnvBasePath = runDirectory + File.separator + PURPLE_DIRECTORY;
        return FittedPurityFile.read(cnvBasePath, sample);
    }

    @NotNull
    static List<PurpleCopyNumber> loadPurpleCopyNumbers(@NotNull final String runDirectory, @NotNull final String sample)
            throws IOException {
        final String cnvBasePath = runDirectory + File.separator + PURPLE_DIRECTORY;
        return PurpleCopyNumberFile.read(cnvBasePath, sample);
    }

    @NotNull
    static List<GeneCopyNumber> loadPurpleGeneCopyNumbers(@NotNull final String runDirectory, @NotNull final String sample)
            throws IOException {
        final String cnvBasePath = runDirectory + File.separator + PURPLE_DIRECTORY;
        final String fileName = GeneCopyNumberFile.generateFilename(cnvBasePath, sample);
        return GeneCopyNumberFile.read(fileName);
    }

    @NotNull
    static Path findStructuralVariantVCF(@NotNull final String runDirectory) throws IOException {
        final Optional<Path> path = Files.walk(Paths.get(runDirectory)).filter(p -> p.toString().endsWith(SV_EXTENSION)).findFirst();
        assert path.isPresent();
        return path.get();
    }

    @NotNull
    static List<SomaticVariant> loadSomaticSNVFile(@NotNull final String sample, @NotNull final String path)
            throws IOException, HartwigException {
        return new SomaticVariantFactory().fromVCFFile(sample, path, SOMATIC_SNV_EXTENSION);
    }

    @NotNull
    static String extractTumorType(@NotNull final CpctEcrfModel cpctEcrfModel, @NotNull final String sample) {
        final String patientId = toPatientId(sample);
        if (patientId == null) {
            LOGGER.warn("Could not resolve patient id from " + sample);
            return Strings.EMPTY;
        }

        final EcrfPatient patient = cpctEcrfModel.findPatientById(patientId);
        if (patient == null) {
            LOGGER.warn("Could not find patient " + patientId + " in CPCT ECRF database!");
            return Strings.EMPTY;
        }

        final List<String> tumorTypesForPatient = readEcrfField(patient, TUMOR_TYPE_ECRF_FIELD);
        if (tumorTypesForPatient == null || tumorTypesForPatient.size() == 0) {
            return Strings.EMPTY;
        }
        // KODU: We should never have more than one tumor type for a single patient.
        assert tumorTypesForPatient.size() == 1;

        if (tumorTypesForPatient.get(0).toLowerCase().startsWith("other")) {
            final List<String> otherTumorTypeForPatient = readEcrfField(patient, TUMOR_TYPE_OTHER_ECRF_FIELD);
            if (otherTumorTypeForPatient == null || otherTumorTypeForPatient.size() == 0) {
                return Strings.EMPTY;
            }
            assert otherTumorTypeForPatient.size() == 1;
            return otherTumorTypeForPatient.get(0);
        } else {
            return tumorTypesForPatient.get(0);
        }
    }

    private static List<String> readEcrfField(@NotNull final EcrfPatient patient, @NotNull final String field) {
        final List<String> fieldValues = patient.fieldValuesByName(field);
        if (fieldValues == null) {
            LOGGER.warn("Could not find field " + field + " in patient " + patient.patientId());
        } else if (fieldValues.size() == 0) {
            LOGGER.warn("No value found for " + field + " in patient " + patient.patientId());
        }
        return fieldValues;
    }

    @Nullable
    private static String toPatientId(@NotNull final String sample) {
        return sample.length() >= 12 ? sample.substring(0, 12) : null;
    }
}
