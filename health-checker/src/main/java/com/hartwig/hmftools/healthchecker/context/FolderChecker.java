package com.hartwig.hmftools.healthchecker.context;

import java.io.File;
import java.io.IOException;
import java.nio.file.DirectoryStream;
import java.nio.file.Files;
import java.nio.file.Path;

import com.hartwig.hmftools.common.exception.EmptyFolderException;
import com.hartwig.hmftools.common.exception.FolderDoesNotExistException;
import com.hartwig.hmftools.common.exception.HealthChecksException;
import com.hartwig.hmftools.common.exception.NotFolderException;

import org.jetbrains.annotations.NotNull;

@FunctionalInterface
public interface FolderChecker {

    @NotNull
    String checkFolder(@NotNull String directory) throws IOException, HealthChecksException;

    @NotNull
    static FolderChecker build() {
        return (directory) -> {
            final File folder = new File(directory);
            checkIfDirectoryExist(folder);
            checkIfIsDirectory(folder);
            checkIfDirectoryIsEmpty(folder);
            return folder.getPath();
        };
    }

    static void checkIfDirectoryExist(@NotNull final File folder) throws FolderDoesNotExistException {
        if (!folder.exists()) {
            throw new FolderDoesNotExistException(folder.toPath().toString());
        }
    }

    static void checkIfIsDirectory(@NotNull final File folder) throws NotFolderException {
        if (!folder.isDirectory()) {
            throw new NotFolderException(folder.toPath().toString());
        }
    }

    static void checkIfDirectoryIsEmpty(@NotNull final File directory) throws IOException, EmptyFolderException {
        final Path path = directory.toPath();
        if (!isDirectoryNotEmpty(path)) {
            throw new EmptyFolderException(path.toString());
        }
    }

    static boolean isDirectoryNotEmpty(@NotNull final Path path) throws IOException {
        try (DirectoryStream<Path> directoryStream = Files.newDirectoryStream(path, getHiddenFilesFilter())) {
            return directoryStream.iterator().hasNext();
        }
    }

    @NotNull
    static DirectoryStream.Filter<Path> getHiddenFilesFilter() {
        return entry -> !Files.isHidden(entry);
    }
}
