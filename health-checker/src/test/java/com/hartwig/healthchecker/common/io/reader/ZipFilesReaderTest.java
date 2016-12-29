package com.hartwig.healthchecker.common.io.reader;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertTrue;

import java.io.FileNotFoundException;
import java.io.IOException;
import java.net.URL;
import java.util.List;

import com.google.common.io.Resources;

import org.junit.Test;

public class ZipFilesReaderTest {

    private static final String ZIP_DIRECTORY = "zipfiles";
    private static final String COMMON_FILE_NAME = "helloworld";
    private static final String HELLO_WORLD_LINE = "helloworld";
    private static final int NUM_HELLO_WORLD_LINES = 2;
    private static final String EMPTY_FILE_NAME = "emptyfile";

    private static final String PATH_WITH_NO_ZIPS = "run";

    @Test
    public void readAllLinesFromZip() throws IOException {
        final URL url = Resources.getResource(ZIP_DIRECTORY);
        final ZipFilesReader zipFileReader = new ZipFilesReader();
        final List<String> lines = zipFileReader.readAllLinesFromZips(url.getPath(), COMMON_FILE_NAME);
        assertEquals(NUM_HELLO_WORLD_LINES, lines.size());
    }

    @Test
    public void readAllLinesFromEmptyFileFromZip() throws IOException {
        final URL url = Resources.getResource(ZIP_DIRECTORY);
        final ZipFilesReader zipFileReader = new ZipFilesReader();
        final List<String> lines = zipFileReader.readAllLinesFromZips(url.getPath(), EMPTY_FILE_NAME);
        assertTrue(lines.isEmpty());
    }

    @Test
    public void readFieldFromZipFiles() throws IOException {
        final URL url = Resources.getResource(ZIP_DIRECTORY);
        final ZipFilesReader zipFileReader = new ZipFilesReader();
        final List<String> lines = zipFileReader.readFieldFromZipFiles(url.getPath(), COMMON_FILE_NAME,
                HELLO_WORLD_LINE);
        assertEquals(NUM_HELLO_WORLD_LINES, lines.size());
    }

    @Test(expected = FileNotFoundException.class)
    public void readFieldFromZipFilesEmpty() throws IOException {
        final URL url = Resources.getResource(PATH_WITH_NO_ZIPS);
        final ZipFilesReader zipFileReader = new ZipFilesReader();
        zipFileReader.readFieldFromZipFiles(url.getPath(), COMMON_FILE_NAME, HELLO_WORLD_LINE);
    }
}
