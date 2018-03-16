package com.hartwig.hmftools.common.io.path;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.Optional;
import java.util.stream.Stream;

import org.jetbrains.annotations.NotNull;

@FunctionalInterface
public interface PathExtensionFinder {

    String FILE_NOT_FOUND_MSG = "File %s not found in path %s";

    @NotNull
    Path findPath(@NotNull String path, @NotNull String extension) throws FileNotFoundException;

    @NotNull
    static PathExtensionFinder build() {
        return (path, extension) -> {
            final Optional<Path> searchedFile = getPath(path, extension);
            if (!searchedFile.isPresent()) {
                throw new FileNotFoundException(String.format(FILE_NOT_FOUND_MSG, extension, path));
            }
            return searchedFile.get();
        };
    }

    @NotNull
    static Optional<Path> getPath(@NotNull final String path, @NotNull final String extension) {
        Stream<Path> paths;
        try {
            paths = Files.walk(new File(path).toPath());
        } catch (IOException e) {
            return Optional.empty();
        }
        return paths.filter(filePath -> filePath.getFileName().toString().endsWith(extension)).findFirst();
    }
}
