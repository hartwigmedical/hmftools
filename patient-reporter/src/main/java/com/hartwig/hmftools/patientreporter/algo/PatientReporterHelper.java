package com.hartwig.hmftools.patientreporter.algo;

import java.io.File;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.List;

import com.hartwig.hmftools.common.ecrf.CpctEcrfModel;
import com.hartwig.hmftools.common.ecrf.datamodel.EcrfPatient;
import com.hartwig.hmftools.common.exception.HartwigException;
import com.hartwig.hmftools.common.gene.GeneCopyNumber;
import com.hartwig.hmftools.common.gene.GeneCopyNumberFile;
import com.hartwig.hmftools.common.purple.copynumber.PurpleCopyNumber;
import com.hartwig.hmftools.common.purple.copynumber.PurpleCopyNumberFile;
import com.hartwig.hmftools.common.purple.purity.FittedPurityFile;
import com.hartwig.hmftools.common.purple.purity.PurityContext;
import com.hartwig.hmftools.common.variant.vcf.VCFFileLoader;
import com.hartwig.hmftools.common.variant.vcf.VCFSomaticFile;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.apache.logging.log4j.util.Strings;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

final class PatientReporterHelper {

    private static final Logger LOGGER = LogManager.getLogger(PatientReporterHelper.class);

    private static final String SOMATIC_EXTENSION = "_post_processed.vcf";
    private static final String PURPLE_DIRECTORY = "purple";
    private static final String MANTA_BPI_FILENAME = "_somaticSV_bpi.vcf";
    private static final String MANTA_FILENAME = "somaticSV.vcf.gz";

    private static final String TUMOR_TYPE_ECRF_FIELD = "BASELINE.CARCINOMA.CARCINOMA.PTUMLOC";

    private PatientReporterHelper() {
    }

    @NotNull
    static PurityContext loadPurity(@NotNull final String runDirectory, @NotNull final String sample) throws IOException, HartwigException {
        final String cnvBasePath = runDirectory + File.separator + PURPLE_DIRECTORY;
        return FittedPurityFile.read(cnvBasePath, sample);
    }

    @NotNull
    static List<PurpleCopyNumber> loadPurpleCopyNumbers(@NotNull final String runDirectory, @NotNull final String sample)
            throws IOException, HartwigException {
        final String cnvBasePath = runDirectory + File.separator + PURPLE_DIRECTORY;
        return PurpleCopyNumberFile.read(cnvBasePath, sample);
    }

    @NotNull
    static List<GeneCopyNumber> loadPurpleGeneCopyNumbers(@NotNull final String runDirectory, @NotNull final String sample)
            throws IOException, HartwigException {
        final String cnvBasePath = runDirectory + File.separator + PURPLE_DIRECTORY;
        final String fileName = GeneCopyNumberFile.generateFilename(cnvBasePath, sample);
        return GeneCopyNumberFile.read(fileName);
    }

    @Nullable
    static Path findMantaVCF(@NotNull final String runDirectory) {
        try {
            final Path vanilla_manta =
                    Files.walk(Paths.get(runDirectory)).filter(p -> p.getFileName().endsWith(MANTA_FILENAME)).findFirst().orElse(null);
            return Files.walk(Paths.get(runDirectory))
                    .filter(p -> p.getFileName().endsWith(MANTA_BPI_FILENAME))
                    .findFirst()
                    .orElse(vanilla_manta);
        } catch (final Exception e) {
            return null;
        }
    }

    @NotNull
    static VCFSomaticFile loadVariantFile(@NotNull final String path) throws IOException, HartwigException {
        return VCFFileLoader.loadSomaticVCF(path, SOMATIC_EXTENSION);
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

        final List<String> tumorTypesForPatient = patient.fieldValuesByName(TUMOR_TYPE_ECRF_FIELD);
        if (tumorTypesForPatient == null) {
            LOGGER.warn("Could not find field " + TUMOR_TYPE_ECRF_FIELD + " in patient " + patientId);
            return Strings.EMPTY;
        }

        if (tumorTypesForPatient.size() == 0) {
            LOGGER.warn("No value found for " + TUMOR_TYPE_ECRF_FIELD + " in patient " + patientId);
            return Strings.EMPTY;
        }

        // KODU: We should never have more than one tumor type for a single patient.
        assert tumorTypesForPatient.size() == 1;
        return tumorTypesForPatient.get(0);
    }

    @Nullable
    private static String toPatientId(@NotNull final String sample) {
        return sample.length() >= 12 ? sample.substring(0, 12) : null;
    }
}
