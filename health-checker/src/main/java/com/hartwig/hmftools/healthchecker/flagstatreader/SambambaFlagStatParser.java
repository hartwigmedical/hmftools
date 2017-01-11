package com.hartwig.hmftools.healthchecker.flagstatreader;

import java.io.File;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.ArrayList;
import java.util.List;
import java.util.Optional;

import com.hartwig.hmftools.common.exception.EmptyFileException;

import org.jetbrains.annotations.NotNull;

public class SambambaFlagStatParser implements FlagStatParser {

    private static final String SEPARATOR_REGEX = " ";
    private static final String FLAGSTAT_SUFFIX = ".flagstat";

    private static final int QC_PASSED_COLUMN = 0;
    private static final int QC_FAILED_COLUMN = 2;
    private static final int FIRST_WORD_COLUMN = 3;

    @Override
    @NotNull
    public FlagStatData parse(@NotNull final String flagstatPath, @NotNull final String filter)
                    throws IOException, EmptyFileException {
        final Optional<Path> filePath = Files.walk(new File(flagstatPath).toPath())
                        .filter(path -> path.getFileName().toString().endsWith(FLAGSTAT_SUFFIX)
                                        && path.getFileName().toString().contains(filter))
                        .findFirst();

        if (!filePath.isPresent()) {
            throw new IOException("Could not find input file with filter \"" + filter + "\" on path " + flagstatPath);
        }

        return parseFromFilePath(filePath.get());
    }

    @NotNull
    private static FlagStatData parseFromFilePath(@NotNull Path filePath) throws IOException, EmptyFileException {
        final int[] index = { 0 };
        final List<FlagStats> passedStats = new ArrayList<>();
        final List<FlagStats> failedStats = new ArrayList<>();

        Files.lines(filePath).map(line -> {
            final String qcPassed = line.split(SEPARATOR_REGEX)[QC_PASSED_COLUMN];
            final String qcFailed = line.split(SEPARATOR_REGEX)[QC_FAILED_COLUMN];

            final String firstWord = line.split(SEPARATOR_REGEX)[FIRST_WORD_COLUMN];
            final int firstWordIndex = line.indexOf(firstWord);
            final String checkName = line.substring(firstWordIndex, line.length());

            return new String[] {qcPassed, qcFailed, checkName};
        }).forEach(array -> {
            final double passedValue = Double.parseDouble(array[0]);
            final double failedValue = Double.parseDouble(array[1]);
            final Optional<FlagStatsType> statsTypeOpt = FlagStatsType.getByIndex(index[0]);
            assert statsTypeOpt.isPresent();

            final FlagStats passed = new FlagStats(statsTypeOpt.get(), passedValue);
            final FlagStats failed = new FlagStats(statsTypeOpt.get(), failedValue);
            passedStats.add(passed);
            failedStats.add(failed);
            index[0] += 1;
        });

        if (passedStats.isEmpty() || failedStats.isEmpty()) {
            throw new EmptyFileException(FLAGSTAT_SUFFIX, filePath.toString());
        }

        return new FlagStatData(passedStats, failedStats);
    }
}
