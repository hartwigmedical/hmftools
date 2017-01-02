package com.hartwig.hmftools.common.io.path;

import static org.junit.Assert.assertNotNull;

import java.io.FileNotFoundException;
import java.io.IOException;
import java.net.URL;
import java.nio.file.NoSuchFileException;
import java.nio.file.Path;

import com.google.common.io.Resources;

import org.junit.Test;

public class PathRegexFinderTest {

    private static final String BASE_DIRECTORY = "PathFinderRegex";
    private static final String EXISTING_FILE = "base_file";
    private static final String EXISTING_REGEX = "file_in_dir";
    private static final String NON_EXISTING_REGEX = "this does not exist";

    private static final String NON_EXISTING_DIRECTORY = "bla";
    private static final String EMPTY_DIRECTORY = "empty";

    @Test
    public void findExactFilePath() throws IOException {
        final URL testPath = Resources.getResource(BASE_DIRECTORY);
        final Path path = PathRegexFinder.build().findPath(testPath.getPath(), EXISTING_FILE);
        assertNotNull(path);
    }

    @Test
    public void findPathOnRegex() throws IOException {
        final URL testPath = Resources.getResource(BASE_DIRECTORY);
        final Path path = PathRegexFinder.build().findPath(testPath.getPath(), EXISTING_REGEX);
        assertNotNull(path);
    }

    @Test(expected = FileNotFoundException.class)
    public void findPathWrongRegex() throws IOException {
        final URL testPath = Resources.getResource(BASE_DIRECTORY);
        PathRegexFinder.build().findPath(testPath.getPath(), NON_EXISTING_REGEX);
    }

    @Test(expected = FileNotFoundException.class)
    public void findPathEmpty() throws IOException {
        final URL testPath = Resources.getResource(EMPTY_DIRECTORY);
        PathRegexFinder.build().findPath(testPath.getPath(), EXISTING_REGEX);
    }

    @Test(expected = NoSuchFileException.class)
    public void findPathNonExistingDir() throws IOException {
        PathRegexFinder.build().findPath(NON_EXISTING_DIRECTORY, EXISTING_REGEX);
    }
}
