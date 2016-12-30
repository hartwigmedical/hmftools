package com.hartwig.healthchecker.io.path;

import static org.junit.Assert.assertTrue;

import java.io.FileNotFoundException;
import java.io.IOException;
import java.nio.file.Path;
import java.util.List;

import com.google.common.io.Resources;

import org.junit.Test;

public class PathsExtensionFinderTest {

    private static final String TEST_DIRECTORY = Resources.getResource("common/zipfiles").getPath();

    private static final String EXISTING_EXTENSION = ".zip";
    private static final String NON_EXISTING_EXTENSION = ".bla";

    @Test
    public void findPaths() throws IOException {
        final List<Path> paths = PathsExtensionFinder.build().findPaths(TEST_DIRECTORY, EXISTING_EXTENSION);
        assertTrue(!paths.isEmpty());
    }

    @Test(expected = FileNotFoundException.class)
    public void findPathsFilesNotFound() throws IOException {
        PathsExtensionFinder.build().findPaths(TEST_DIRECTORY, NON_EXISTING_EXTENSION);
    }
}
