package com.hartwig.hmftools.common.copynumber.freec;

import com.hartwig.hmftools.common.exception.EmptyFileException;
import com.hartwig.hmftools.common.exception.HartwigException;
import com.hartwig.hmftools.common.io.path.PathRegexFinder;
import com.hartwig.hmftools.common.io.reader.FileReader;
import com.hartwig.hmftools.common.io.reader.LineReader;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

import java.io.File;
import java.io.IOException;
import java.nio.file.Path;
import java.util.Collections;
import java.util.List;

public final class FreecFileLoader {

    private static final Logger LOGGER = LogManager.getLogger(FreecFileLoader.class);

    // KODU: copynumber data is stored in {run}/copyNumber/{sampleR}_{sampleT}/freec/{sampleT}<>.bam_CNVs
    private static final String COPYNUMBER_BASE_DIRECTORY = "copyNumber";
    private static final String COPYNUMBER_SAMPLE_CONNECTOR = "_";
    private static final String FREEC_ALGO_DIRECTORY = "freec";

    private static final String COPYNUMBER_REGEX = "%s[_.].*(?<!_normal)_CNVs$";

    private static final String NORMAL_RATIO_REGEX = "%s[_.].*_normal_ratio.txt$";
    private static final String TUMOR_RATIO_REGEX = "%s[_.].*(?<!_normal)_ratio.txt$";

    private static final String RATIO_START_OF_HEADER = "Chromosome";

    @NotNull
    public static String getFreecBasePath(@NotNull final String runDirectory, @NotNull final String refSample,
            @Nullable final String tumorSample) {
        final String baseDir = runDirectory + File.separator + COPYNUMBER_BASE_DIRECTORY + File.separator;
        final String setDir = tumorSample != null ? refSample + COPYNUMBER_SAMPLE_CONNECTOR + tumorSample : refSample;

        return baseDir + setDir + File.separator + FREEC_ALGO_DIRECTORY;
    }

    @NotNull
    static List<String> copyNumberLines(@NotNull final String basePath, @NotNull final String sample)
            throws EmptyFileException, IOException {
        final Path copynumberPath = PathRegexFinder.build().findPath(basePath,
                String.format(COPYNUMBER_REGEX, sample));
        try {
            LOGGER.debug("Loading CNV data from " + copynumberPath.toString());
            return FileReader.build().readLines(copynumberPath);
        } catch (EmptyFileException e) {
            // if the CNV is empty (but exists) and the ratio file exists, there is no problem (just no CNVs found)
            final Path ratioPath = PathRegexFinder.build().findPath(basePath,
                    String.format(TUMOR_RATIO_REGEX, sample));
            FileReader.build().readLines(ratioPath);
            return Collections.emptyList();
        }
    }

    @NotNull
    static List<String> normalRatioLines(@NotNull final String basePath, @NotNull final String sample)
            throws HartwigException, IOException {
        return ratioLines(basePath, sample, NORMAL_RATIO_REGEX);
    }

    @NotNull
    static List<String> tumorRatioLines(@NotNull final String basePath, @NotNull final String sample)
            throws HartwigException, IOException {
        return ratioLines(basePath, sample, TUMOR_RATIO_REGEX);
    }

    @NotNull
    private static List<String> ratioLines(@NotNull final String basePath, @NotNull final String sample,
            @NotNull final String regex) throws HartwigException, IOException {
        final Path path = PathRegexFinder.build().findPath(basePath, String.format(regex, sample));

        LOGGER.debug("Loading ratio data from " + path.toString());
        return LineReader.build().readLines(path, x -> !x.startsWith(RATIO_START_OF_HEADER));
    }
}
