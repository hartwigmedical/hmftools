package com.hartwig.hmftools.common.io.path;

import static org.junit.Assert.assertTrue;

import java.io.FileNotFoundException;
import java.nio.file.Path;
import java.util.List;

import com.google.common.io.Resources;

import org.junit.Test;

public class PathsExtensionFinderTest {

    private static final String TEST_DIRECTORY = Resources.getResource("io").getPath();

    private static final String EXISTING_EXTENSION = ".file";
    private static final String NON_EXISTING_EXTENSION = ".bla";

    @Test
    public void findCorrectPaths() throws FileNotFoundException {
        final List<Path> paths = PathsExtensionFinder.build().findPaths(TEST_DIRECTORY, EXISTING_EXTENSION);
        assertTrue(!paths.isEmpty());
    }

    @Test(expected = FileNotFoundException.class)
    public void throwErrorWhenExtensionDoesNotExist() throws FileNotFoundException {
        PathsExtensionFinder.build().findPaths(TEST_DIRECTORY, NON_EXISTING_EXTENSION);
    }
}
