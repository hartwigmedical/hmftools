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
import com.hartwig.hmftools.common.exception.EmptyFileException;
import com.hartwig.hmftools.common.exception.HartwigException;
import com.hartwig.hmftools.common.variant.vcf.VCFFileLoader;
import com.hartwig.hmftools.common.variant.vcf.VCFSomaticFile;

import org.jetbrains.annotations.NotNull;

final class PatientReporterHelper {

    private static final String SOMATIC_EXTENSION = "_melted.vcf";
    private static final String COPYNUMBER_DIRECTORY = "copyNumber";
    private static final String COPYNUMBER_EXTENSION = ".bam_CNVs";
    private static final String FREEC_DIRECTORY = "freec";

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
