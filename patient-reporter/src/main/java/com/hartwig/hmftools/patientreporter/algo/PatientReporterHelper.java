package com.hartwig.hmftools.patientreporter.algo;

import java.io.File;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.List;
import java.util.Optional;
import java.util.stream.Collectors;

import com.hartwig.hmftools.common.ecrf.projections.PatientCancerTypes;
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
    private static final String CIRCOS_PLOT_DIRECTORY = "plot";
    private static final String CIRCOS_PLOT_EXTENSION = ".circos.png";

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
    static String findCircosPlotPath(@NotNull final String runDirectory, @NotNull final String sample) {
        return runDirectory + File.separator + PURPLE_DIRECTORY + File.separator + CIRCOS_PLOT_DIRECTORY + File.separator + sample
                + CIRCOS_PLOT_EXTENSION;
    }

    @NotNull
    static List<SomaticVariant> loadPassedSomaticVariants(@NotNull final String sample, @NotNull final String path) throws IOException {
        return SomaticVariantFactory.passOnlyInstance().fromVCFFile(sample, path, SOMATIC_SNV_EXTENSION);
    }

    @NotNull
    static String extractCancerType(@NotNull final List<PatientCancerTypes> patientsCancerTypes, @NotNull final String sample) {
        final String patientId = toPatientId(sample);
        if (patientId == null) {
            LOGGER.warn("Could not resolve patient id from " + sample);
            return Strings.EMPTY;
        }
        final List<PatientCancerTypes> matchingIdCancerTypes = patientsCancerTypes.stream()
                .filter(patientCancerTypes -> patientCancerTypes.cpctId().equals(patientId))
                .collect(Collectors.toList());

        // KODU: We should never have more than one cancer type for a single patient.
        assert matchingIdCancerTypes.size() < 2;

        if (matchingIdCancerTypes.size() == 1) {
            return matchingIdCancerTypes.get(0).cancerType();
        } else {
            LOGGER.warn("Could not find patient " + patientId + " in CPCT ECRF database!");
            return Strings.EMPTY;
        }
    }

    @Nullable
    private static String toPatientId(@NotNull final String sample) {
        return sample.length() >= 12 ? sample.substring(0, 12) : null;
    }
}
