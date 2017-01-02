package com.hartwig.hmftools.common.io.reader;

import static java.util.stream.Collectors.toList;

import java.io.IOException;
import java.nio.file.Path;
import java.util.ArrayList;
import java.util.Collection;
import java.util.List;

import com.hartwig.hmftools.common.exception.HealthChecksException;
import com.hartwig.hmftools.common.io.path.PathsExtensionFinder;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

public class ZipFilesReader {

    private static final String ERROR_MSG =
            "Error occurred when reading file %s. " + "Will return empty stream. Error -> %s";

    private static final Logger LOGGER = LogManager.getLogger(ZipFilesReader.class);
    private static final String ZIP_EXTENSION = ".zip";

    @NotNull
    public List<String> readAllLinesFromZips(@NotNull final String path, @NotNull final String fileName)
            throws IOException {
        final List<Path> zipPaths = PathsExtensionFinder.build().findPaths(path, ZIP_EXTENSION);
        return zipPaths.stream().map(zipPath -> readFileFromZip(zipPath.toString(), fileName)).flatMap(
                Collection::stream).collect(toList());
    }

    @NotNull
    public List<String> readFieldFromZipFiles(@NotNull final String path, @NotNull final String fileName,
            @NotNull final String filter) throws IOException {
        final List<Path> zipPaths = PathsExtensionFinder.build().findPaths(path, ZIP_EXTENSION);
        return zipPaths.stream().map(zipPath -> searchForLineInZip(zipPath, fileName, filter)).filter(
                line -> line != null).collect(toList());
    }

    @NotNull
    private static List<String> readFileFromZip(@NotNull final String path, @NotNull final String fileName) {
        final List<String> fileLines = new ArrayList<>();
        try {
            fileLines.addAll(FileInZipsReader.build().readLines(path, fileName));
        } catch (IOException | HealthChecksException e) {
            LOGGER.error(String.format(ERROR_MSG, path, e.getMessage()));
        }
        return fileLines;
    }

    @Nullable
    private static String searchForLineInZip(@NotNull final Path path, @NotNull final String fileName,
            @NotNull final String filter) {
        String searchedLine = null;
        try {
            searchedLine = LineInZipsReader.build().readLines(path.toString(), fileName, filter);
        } catch (IOException | HealthChecksException e) {
            LOGGER.error(String.format(ERROR_MSG, path, e.getMessage()));
        }
        return searchedLine;
    }
}
