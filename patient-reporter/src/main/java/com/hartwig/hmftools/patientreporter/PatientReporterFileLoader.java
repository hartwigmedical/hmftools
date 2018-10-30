package com.hartwig.hmftools.patientreporter;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.List;
import java.util.Optional;
import java.util.stream.Collectors;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.io.path.PathExtensionFinder;
import com.hartwig.hmftools.common.purple.copynumber.PurpleCopyNumber;
import com.hartwig.hmftools.common.purple.copynumber.PurpleCopyNumberFile;
import com.hartwig.hmftools.common.purple.gene.GeneCopyNumber;
import com.hartwig.hmftools.common.purple.gene.GeneCopyNumberFile;
import com.hartwig.hmftools.common.purple.purity.FittedPurityFile;
import com.hartwig.hmftools.common.purple.purity.PurityContext;
import com.hartwig.hmftools.common.variant.SomaticVariant;
import com.hartwig.hmftools.common.variant.SomaticVariantFactory;
import com.hartwig.hmftools.common.variant.enrich.SomaticEnrichment;
import com.hartwig.hmftools.patientreporter.chord.ChordAnalysis;
import com.hartwig.hmftools.patientreporter.chord.ChordFile;
import com.hartwig.hmftools.patientreporter.germline.BachelorFile;
import com.hartwig.hmftools.patientreporter.germline.GermlineVariant;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

import htsjdk.variant.variantcontext.filter.PassingVariantFilter;

public final class PatientReporterFileLoader {

    private static final Logger LOGGER = LogManager.getLogger(PatientReporterFileLoader.class);

    private static final String PURPLE_DIRECTORY = "purple";
    private static final String CIRCOS_PLOT_DIRECTORY = "plot";
    private static final String CIRCOS_PLOT_EXTENSION = ".circos.png";
    private static final String SV_EXTENSION_V3 = "_somaticSV_bpi.vcf";
    private static final String SV_EXTENSION_V4 = "_somaticSV_bpi.vcf.gz";
    private static final String SOMATIC_VCF_EXTENSION_V3 = "_post_processed_v2.2.vcf.gz";
    private static final String SOMATIC_VCF_EXTENSION_V4 = "_post_processed.vcf.gz";
    private static final String BACHELOR_DIRECTORY = "bachelor";
    private static final String CHORD_DIRECTORY = "chord";

    private PatientReporterFileLoader() {
    }

    @NotNull
    static PurityContext loadPurity(@NotNull String runDirectory, @NotNull String sample) throws IOException {
        final String cnvBasePath = runDirectory + File.separator + PURPLE_DIRECTORY;
        return FittedPurityFile.read(cnvBasePath, sample);
    }

    @NotNull
    static List<PurpleCopyNumber> loadPurpleCopyNumbers(@NotNull String runDirectory, @NotNull String sample) throws IOException {
        final String cnvBasePath = runDirectory + File.separator + PURPLE_DIRECTORY;
        return PurpleCopyNumberFile.read(PurpleCopyNumberFile.generateFilename(cnvBasePath, sample));
    }

    @NotNull
    static List<GeneCopyNumber> loadPurpleGeneCopyNumbers(@NotNull String runDirectory, @NotNull String sample) throws IOException {
        final String cnvBasePath = runDirectory + File.separator + PURPLE_DIRECTORY;
        final String fileName = GeneCopyNumberFile.generateFilename(cnvBasePath, sample);
        return GeneCopyNumberFile.read(fileName);
    }

    @NotNull
    static String findCircosPlotPath(@NotNull String runDirectory, @NotNull String sample) {
        return runDirectory + File.separator + PURPLE_DIRECTORY + File.separator + CIRCOS_PLOT_DIRECTORY + File.separator + sample
                + CIRCOS_PLOT_EXTENSION;
    }

    @NotNull
    static Path findStructuralVariantVCF(@NotNull String runDirectory) throws IOException {
        // TODO (KODU): Clean up once pipeline v3 no longer exists
        Optional<Path> path = Files.walk(Paths.get(runDirectory)).filter(p -> p.toString().endsWith(SV_EXTENSION_V3)).findFirst();
        if (!path.isPresent()) {
            path = Files.walk(Paths.get(runDirectory)).filter(p -> p.toString().endsWith(SV_EXTENSION_V4)).findFirst();
        }
        assert path.isPresent();
        return path.get();
    }

    @NotNull
    static List<SomaticVariant> loadPassedSomaticVariants(@NotNull String runDirectory, @NotNull String sample,
            @NotNull SomaticEnrichment somaticEnrichment) throws IOException {
        // TODO (KODU): Clean up once pipeline v3 no longer exists
        Path vcfPath;
        try {
            vcfPath = PathExtensionFinder.build().findPath(runDirectory, SOMATIC_VCF_EXTENSION_V3);
        } catch (FileNotFoundException exception) {
            vcfPath = PathExtensionFinder.build().findPath(runDirectory, SOMATIC_VCF_EXTENSION_V4);
        }

        return SomaticVariantFactory.filteredInstanceWithEnrichment(new PassingVariantFilter(), somaticEnrichment)
                .fromVCFFile(sample, vcfPath.toString());
    }

    @Nullable
    static List<GermlineVariant> loadPassedGermlineVariants(@NotNull String runDirectory, @NotNull String sample) throws IOException {
        String bachelorDirectory = runDirectory + File.separator + BACHELOR_DIRECTORY;
        if (!BachelorFile.hasBachelorRun(bachelorDirectory, sample)) {
            return null;
        } else {
            String bachelorFile = BachelorFile.findBachelorFilePath(bachelorDirectory, sample);
            List<GermlineVariant> germlineVariants =
                    bachelorFile != null ? BachelorFile.loadBachelorFile(bachelorFile) : Lists.newArrayList();

            return germlineVariants.stream().filter(GermlineVariant::passFilter).collect(Collectors.toList());
        }
    }

    @NotNull
    static ChordAnalysis loadChordFile(@NotNull String runDirectory, @NotNull String sample) throws IOException {
        final String chordDirectory = runDirectory + File.separator + CHORD_DIRECTORY;
        return ChordFile.loadChordFile(chordDirectory, sample);
    }
}
