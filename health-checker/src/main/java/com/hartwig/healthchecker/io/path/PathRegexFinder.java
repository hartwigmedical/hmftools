package com.hartwig.healthchecker.io.path;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.Optional;
import java.util.stream.Stream;

import org.jetbrains.annotations.NotNull;

@FunctionalInterface
public interface PathRegexFinder {

    String FILE_S_NOT_FOUND_MSG = "No match for %s found in path %s";

    @NotNull
    Path findPath(@NotNull String path, @NotNull String regex) throws IOException;

    @NotNull
    static PathRegexFinder build() {
        return (path, regex) -> {
            final Optional<Path> searchedFile = getPath(path, regex);
            if (!searchedFile.isPresent()) {
                throw new FileNotFoundException(String.format(FILE_S_NOT_FOUND_MSG, regex, path));
            }
            return searchedFile.get();
        };
    }

    @NotNull
    static Optional<Path> getPath(@NotNull final String path, @NotNull final String regex) throws IOException {
        try (Stream<Path> paths = Files.walk(new File(path).toPath())) {
            return paths.filter(filePath -> filePath.getFileName().toString().matches(regex)).findFirst();
        }
    }
}
