package com.hartwig.hmftools.patientreporter;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.List;
import java.util.stream.Collectors;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.copynumber.CopyNumber;
import com.hartwig.hmftools.common.copynumber.cnv.CNVFileLoader;
import com.hartwig.hmftools.common.ecrf.CpctEcrfModel;
import com.hartwig.hmftools.common.ecrf.datamodel.EcrfField;
import com.hartwig.hmftools.common.ecrf.datamodel.EcrfPatient;
import com.hartwig.hmftools.common.exception.EmptyFileException;
import com.hartwig.hmftools.common.exception.HartwigException;
import com.hartwig.hmftools.common.variant.vcf.VCFFileLoader;
import com.hartwig.hmftools.common.variant.vcf.VCFSomaticFile;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.apache.logging.log4j.util.Strings;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

final class PatientReporterHelper {

    private static final Logger LOGGER = LogManager.getLogger(PatientReporterHelper.class);

    private static final String SOMATIC_EXTENSION = "_melted.vcf";
    private static final String COPYNUMBER_DIRECTORY = "copyNumber";
    private static final String COPYNUMBER_EXTENSION = ".bam_CNVs";
    private static final String FREEC_DIRECTORY = "freec";

    private static final String TUMOR_TYPE_ECRF_FIELD = "BASELINE.CARCINOMA.CARCINOMA.PTUMLOC";

    private PatientReporterHelper() {
    }

    @NotNull
    static VCFSomaticFile loadVariantFile(@NotNull final String path) throws IOException, HartwigException {
        return VCFFileLoader.loadSomaticVCF(path, SOMATIC_EXTENSION);
    }

    @NotNull
    static List<CopyNumber> loadCNVFile(@NotNull final String runDirectory, @NotNull final String sample)
            throws IOException, HartwigException {
        final String cnvBasePath = guessCNVBasePath(runDirectory, sample) + File.separator + FREEC_DIRECTORY;

        try {
            return CNVFileLoader.loadCNV(cnvBasePath, sample, COPYNUMBER_EXTENSION);
        } catch (EmptyFileException e) {
            // KODU: It could be that the sample simply does not have any amplifications...
            return Lists.newArrayList();
        }
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

        final EcrfField tumorTypeField = cpctEcrfModel.findFieldById(TUMOR_TYPE_ECRF_FIELD);
        if (tumorTypeField == null) {
            LOGGER.warn("Could not find field " + TUMOR_TYPE_ECRF_FIELD + " in CPCT ECRF database!");
            return Strings.EMPTY;
        }

        final List<String> valueForPatient = patient.fieldValues(tumorTypeField);
        if (valueForPatient == null) {
            LOGGER.warn("Could not find field " + TUMOR_TYPE_ECRF_FIELD + " in patient " + patientId);
            return Strings.EMPTY;
        }

        if (valueForPatient.size() != 1) {
            LOGGER.warn("No value found for " + TUMOR_TYPE_ECRF_FIELD + " in patient " + patientId);
            return Strings.EMPTY;
        }

        return valueForPatient.get(0);
    }

    @Nullable
    private static String toPatientId(@NotNull final String sample) {
        return sample.length() >= 12 ? sample.substring(0, 12) : null;
    }

    @NotNull
    private static String guessCNVBasePath(@NotNull final String runDirectory, @NotNull final String sample)
            throws IOException {
        final String basePath = runDirectory + File.separator + COPYNUMBER_DIRECTORY;

        for (final Path path : Files.list(new File(basePath).toPath()).collect(Collectors.toList())) {
            if (path.toFile().isDirectory() && path.getFileName().toFile().getName().contains(sample)) {
                return path.toString();
            }
        }

        throw new FileNotFoundException(
                "Could not determine CNV location in " + runDirectory + " using sample " + sample);
    }
}
