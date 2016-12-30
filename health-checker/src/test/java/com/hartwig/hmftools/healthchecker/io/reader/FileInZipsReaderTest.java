package com.hartwig.hmftools.healthchecker.io.reader;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertNotNull;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.net.URL;
import java.util.List;

import com.google.common.io.Resources;
import com.hartwig.hmftools.healthchecker.exception.EmptyFileException;
import com.hartwig.hmftools.healthchecker.exception.HealthChecksException;

import org.junit.Test;

public class FileInZipsReaderTest {

    private static final String ZIP_DIRECTORY = "zipfiles";

    private static final String FILE = "helloworld";
    private static final String FILE_ZIPPED = FILE + ".zip";
    private static final int FILE_EXPECTED_NUM_LINES = 1;

    private static final String EMPTY_FILE = "emptyfile";
    private static final String EMPTY_FILE_ZIPPED = EMPTY_FILE + ".zip";

    @Test
    public void readLines() throws IOException, HealthChecksException {
        final URL testPath = Resources.getResource(ZIP_DIRECTORY + File.separator + FILE_ZIPPED);
        final List<String> readLines = FileInZipsReader.build().readLines(testPath.getPath(), FILE);
        assertNotNull(readLines);
        assertEquals(FILE_EXPECTED_NUM_LINES, readLines.size());
    }

    @Test(expected = EmptyFileException.class)
    public void readLinesEmptyFiles() throws IOException, HealthChecksException {
        final URL testPath = Resources.getResource(ZIP_DIRECTORY + File.separator + EMPTY_FILE_ZIPPED);
        FileInZipsReader.build().readLines(testPath.getPath(), EMPTY_FILE);
    }

    @Test(expected = FileNotFoundException.class)
    public void readLinesFileNotInZip() throws IOException, HealthChecksException {
        final URL testPath = Resources.getResource(ZIP_DIRECTORY + File.separator + FILE_ZIPPED);
        FileInZipsReader.build().readLines(testPath.getPath(), "bla");
    }
}
