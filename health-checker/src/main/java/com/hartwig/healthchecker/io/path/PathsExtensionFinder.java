package com.hartwig.healthchecker.io.path;

import static java.util.stream.Collectors.toCollection;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.ArrayList;
import java.util.List;
import java.util.stream.Stream;

import org.jetbrains.annotations.NotNull;

@FunctionalInterface
public interface PathsExtensionFinder {

    String FILE_S_NOT_FOUND_MSG = "File %s not found in path %s";

    @NotNull
    List<Path> findPaths(@NotNull String path, @NotNull String extension) throws IOException;

    @NotNull
    static PathsExtensionFinder build() {
        return (path, extension) -> {
            final List<Path> searchedFile = getPath(path, extension);
            if (searchedFile.isEmpty()) {
                throw new FileNotFoundException(String.format(FILE_S_NOT_FOUND_MSG, extension, path));
            }
            return searchedFile;
        };
    }

    @NotNull
    static List<Path> getPath(@NotNull final String path, @NotNull final String extension) throws IOException {
        try (Stream<Path> paths = Files.walk(new File(path).toPath())) {
            return paths.filter(filePath -> filePath.getFileName().toString().endsWith(extension)).sorted()
                            .collect(toCollection(ArrayList<Path>::new));
        }
    }
}
