package com.hartwig.hmftools.healthchecker.io.reader;

import static org.junit.Assert.assertTrue;

import java.io.File;
import java.io.IOException;
import java.net.URL;

import com.google.common.io.Resources;
import com.hartwig.hmftools.healthchecker.exception.HealthChecksException;
import com.hartwig.hmftools.healthchecker.exception.LineNotFoundException;

import org.junit.Test;

public class LineInZipsReaderTest {

    private static final String ZIP_DIRECTORY = "zipfiles";
    private static final String FILE = "helloworld";
    private static final String FILE_CONTENT = "helloworld";
    private static final String ZIP_FILE = FILE + ".zip";

    private static final String EMPTY_FILE = "emptyfile";
    private static final String EMPTY_ZIP = EMPTY_FILE + ".zip";

    @Test
    public void readLine() throws IOException, HealthChecksException {
        final URL testPath = Resources.getResource(ZIP_DIRECTORY + File.separator + ZIP_FILE);
        final String line = LineInZipsReader.build().readLines(testPath.getPath(), FILE, FILE_CONTENT);
        assertTrue(line.startsWith(FILE_CONTENT));
        assertTrue(line.endsWith(FILE_CONTENT));
    }

    @Test(expected = LineNotFoundException.class)
    public void readLinesEmptyFiles() throws IOException, HealthChecksException {
        final URL testPath = Resources.getResource(ZIP_DIRECTORY + File.separator + EMPTY_ZIP);
        LineInZipsReader.build().readLines(testPath.getPath(), EMPTY_FILE, "try to find me");
    }

    @Test(expected = LineNotFoundException.class)
    public void readLineNotInFile() throws IOException, HealthChecksException {
        final URL testPath = Resources.getResource(ZIP_DIRECTORY + File.separator + ZIP_FILE);
        LineInZipsReader.build().readLines(testPath.getPath(), FILE, "try to find me");
    }
}
