package com.hartwig.hmftools.common.copynumber.cnv;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.copynumber.CopyNumber;
import com.hartwig.hmftools.common.copynumber.CopyNumberFactory;
import com.hartwig.hmftools.common.exception.EmptyFileException;
import com.hartwig.hmftools.common.exception.HartwigException;
import com.hartwig.hmftools.common.io.path.PathRegexFinder;
import com.hartwig.hmftools.common.io.reader.FileReader;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.jetbrains.annotations.NotNull;

import java.io.IOException;
import java.nio.file.Path;
import java.util.Collections;
import java.util.List;

public final class CNVFileLoader {

    private static final Logger LOGGER = LogManager.getLogger(CNVFileLoader.class);

    private static final String COPYNUMBER_REGEX = "%s[_.].*(?<!_normal)_CNVs$";
    private static final String COPYNUMBER_RATIO_REGEX = "%s[_.].*(?<!_normal)_ratio.txt$";

    @NotNull
    public static List<CopyNumber> loadCNV(@NotNull final String basePath, @NotNull final String sample) throws IOException, HartwigException {
        final List<CopyNumber> copyNumbers = Lists.newArrayList();
        for (final String line : copyNumberLines(basePath, sample)) {
            copyNumbers.add(CopyNumberFactory.fromCNVLine(line));
        }
        return copyNumbers;
    }

    @NotNull
    public static String copyNumberPath(@NotNull final String basePath, @NotNull final String sample) throws IOException {
        return PathRegexFinder.build()
                .findPath(basePath, String.format(COPYNUMBER_REGEX, sample)).toString();
    }

    @NotNull
    private static List<String> copyNumberLines(@NotNull final String basePath, @NotNull final String sample)
            throws EmptyFileException, IOException {
        final Path copynumberPath = PathRegexFinder.build()
                                                   .findPath(basePath, String.format(COPYNUMBER_REGEX, sample));
        try {
            LOGGER.debug("Loading CNV data from " + copynumberPath.toString());
            return FileReader.build().readLines(copynumberPath);
        } catch (EmptyFileException e) {
            // if the CNV is empty (but exists) and the ratio file exists, there is no problem (just no CNVs found)
            // JOBA: Not sure how this works... ratio file contains a header line and CopyNumber is a different field?!?
            final Path ratioPath = PathRegexFinder.build()
                                                  .findPath(basePath, String.format(COPYNUMBER_RATIO_REGEX, sample));
            FileReader.build().readLines(ratioPath);
            return Collections.emptyList();
        }
    }
}
