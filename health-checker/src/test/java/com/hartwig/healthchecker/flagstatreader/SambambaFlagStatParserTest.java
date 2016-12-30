package com.hartwig.healthchecker.flagstatreader;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertNotNull;

import java.io.IOException;
import java.util.List;
import java.util.Optional;

import com.google.common.io.Resources;
import com.hartwig.healthchecker.exception.EmptyFileException;

import org.jetbrains.annotations.NotNull;
import org.junit.Test;

public class SambambaFlagStatParserTest {

    private static final String FLAGSTAT_EXAMPLE_PATH = Resources.getResource("common/flagstat").getPath();
    private static final String FLAGSTAT_EXAMPLE_FILE_PATTERN = ".example";
    private static final String FLAGSTAT_EMPTY_FILE_PATTERN = ".empty";
    private static final String FLAGSTAT_DOES_NOT_EXIST_PATTERN = "this does not exist";

    @Test
    public void parseExampleFile() throws IOException, EmptyFileException {
        final FlagStatParser parser = new SambambaFlagStatParser();
        final FlagStatData flagStatData = parser.parse(FLAGSTAT_EXAMPLE_PATH, FLAGSTAT_EXAMPLE_FILE_PATTERN);

        assertFlagStatData(flagStatData.getPassedStats(), 13, 1D);
        assertFlagStatData(flagStatData.getFailedStats(), 13, 20D);
    }

    @Test(expected = IOException.class)
    public void parseFileNotFound() throws IOException, EmptyFileException {
        final FlagStatParser parser = new SambambaFlagStatParser();
        parser.parse(FLAGSTAT_EXAMPLE_PATH, FLAGSTAT_DOES_NOT_EXIST_PATTERN);
    }

    @Test(expected = EmptyFileException.class)
    public void parseEmptyFile() throws IOException, EmptyFileException {
        final FlagStatParser parser = new SambambaFlagStatParser();
        parser.parse(FLAGSTAT_EXAMPLE_PATH, FLAGSTAT_EMPTY_FILE_PATTERN);
    }

    private static void assertFlagStatData(@NotNull final List<FlagStats> flagStat, final int expectedSize,
            final double expectedTotalIndex) {
        assertEquals(expectedSize, flagStat.size());

        final Optional<FlagStats> passedFlagStat = flagStat.stream().filter(
                flagStats -> flagStats.getFlagStatsType() == FlagStatsType.TOTAL_INDEX).findFirst();

        assert passedFlagStat.isPresent();

        assertNotNull(passedFlagStat.get());
        assertEquals(expectedTotalIndex, passedFlagStat.get().getValue(), 0.0d);
    }
}
