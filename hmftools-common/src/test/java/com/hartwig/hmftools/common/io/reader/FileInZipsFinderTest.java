package com.hartwig.hmftools.common.io.reader;

import static org.junit.Assert.assertNotNull;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.net.URL;
import java.util.List;
import java.util.zip.ZipEntry;
import java.util.zip.ZipException;
import java.util.zip.ZipFile;

import com.google.common.io.Resources;
import com.hartwig.hmftools.common.exception.HealthChecksException;

import org.junit.Test;

public class FileInZipsFinderTest {

    private static final String EMPTY_FILE = "emptyfile";
    private static final String ZIP_WITH_EMPTY_FILE = "emptyfile.zip";
    private static final String ZIP_WITHOUT_FILE = "emptyarchive.zip";
    private static final String ZIP_DIRECTORY = "zipfiles";

    @Test
    public void findFile() throws IOException, HealthChecksException {
        final URL testPath = Resources.getResource(ZIP_DIRECTORY + File.separator + ZIP_WITH_EMPTY_FILE);
        final ZipFile zipFile = new ZipFile(testPath.getPath());
        final List<? extends ZipEntry> zipEntry = FileInZipsFinder.build().findFileInZip(zipFile, EMPTY_FILE);
        assertNotNull(zipEntry);
    }

    @Test(expected = ZipException.class)
    public void readZipNotFound() throws IOException, HealthChecksException {
        final URL testPath = Resources.getResource(ZIP_DIRECTORY + File.separator + ZIP_WITHOUT_FILE);
        final ZipFile zipFile = new ZipFile(testPath.getPath());
        FileInZipsFinder.build().findFileInZip(zipFile, "bla");
    }

    @Test(expected = FileNotFoundException.class)
    public void readFileNotFound() throws IOException, HealthChecksException {
        final URL testPath = Resources.getResource(ZIP_DIRECTORY + File.separator + ZIP_WITH_EMPTY_FILE);
        final ZipFile zipFile = new ZipFile(testPath.getPath());
        FileInZipsFinder.build().findFileInZip(zipFile, "bla");
    }
}
