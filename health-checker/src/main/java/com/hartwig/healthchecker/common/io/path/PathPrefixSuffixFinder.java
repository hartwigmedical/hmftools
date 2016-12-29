package com.hartwig.healthchecker.common.io.path;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.Optional;
import java.util.stream.Stream;

import org.jetbrains.annotations.NotNull;

@FunctionalInterface
public interface PathPrefixSuffixFinder {

    String FILE_S_NOT_FOUND_MSG = "File %s not found in path %s";

    @NotNull
    Path findPath(@NotNull String path, @NotNull String prefix, @NotNull String suffix)
            throws IOException;

    @NotNull
    static PathPrefixSuffixFinder build() {
        return (path, prefix, suffix) -> {
            final Optional<Path> fileFound = getPath(path, prefix, suffix);
            if (!fileFound.isPresent()) {
                throw new FileNotFoundException(String.format(FILE_S_NOT_FOUND_MSG, suffix, path));
            }
            return fileFound.get();
        };
    }

    @NotNull
    static Optional<Path> getPath(@NotNull final String path, @NotNull final String prefix,
            @NotNull final String suffix) throws IOException {
        try (Stream<Path> paths = Files.walk(new File(path).toPath())) {
            return paths.filter(
                    filePath -> filePath.getFileName().toString().startsWith(prefix)
                            && filePath.getFileName().toString().endsWith(suffix)
                            && filePath.toString().contains(path + File.separator + prefix))
                    .findFirst();
        }
    }
}
