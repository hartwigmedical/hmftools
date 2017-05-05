package com.hartwig.hmftools.common.ratio.txt;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.exception.HartwigException;
import com.hartwig.hmftools.common.io.path.PathRegexFinder;
import com.hartwig.hmftools.common.io.reader.LineReader;
import com.hartwig.hmftools.common.ratio.Ratio;
import com.hartwig.hmftools.common.ratio.RatioFactory;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.jetbrains.annotations.NotNull;

import java.io.IOException;
import java.nio.file.Path;
import java.util.Collections;
import java.util.List;

public final class RatioFileLoader {

    private static final Logger LOGGER = LogManager.getLogger(RatioFileLoader.class);
    private static final String START_OF_HEADER = "Chromosome";

    private static final String NORMAL_RATIO_REGEX = "%s[_.].*_normal_ratio.txt$";
    private static final String TUMOR_RATIO_REGEX = "%s[_.].*(?<!_normal)_ratio.txt$";

    @NotNull
    public static List<Ratio> loadNormalRatios(@NotNull final String basePath, @NotNull final String sample) throws IOException, HartwigException {
        return loadRatios(basePath, sample, NORMAL_RATIO_REGEX);
    }

    @NotNull
    public static List<Ratio> loadTumorRatios(@NotNull final String basePath, @NotNull final String sample) throws IOException, HartwigException {
        return loadRatios(basePath, sample, TUMOR_RATIO_REGEX);
    }

    @NotNull
    private static List<Ratio> loadRatios(@NotNull final String basePath, @NotNull final String sample, @NotNull final String regex) throws IOException, HartwigException {
        final List<Ratio> results = Lists.newArrayList();
        for (final String line : ratioLines(basePath, sample, regex)) {
            Ratio ratio = RatioFactory.fromRatioLine(line);
            if (ratio.ratio() > -1) {
                results.add(ratio);
            }
        }

        Collections.sort(results);
        return results;
    }

    @NotNull
    private static List<String> ratioLines(@NotNull final String basePath, @NotNull final String sample, @NotNull final String regex)
            throws HartwigException, IOException {
        final Path path = PathRegexFinder.build().findPath(basePath, String.format(regex, sample));

        LOGGER.debug("Loading ratio data from " + path.toString());
        return LineReader.build().readLines(path, s -> !s.startsWith(START_OF_HEADER));
    }
}
