package com.hartwig.hmftools.healthchecker.context;

import java.io.File;
import java.io.IOException;
import java.nio.file.DirectoryStream;
import java.nio.file.Files;
import java.nio.file.Path;

import com.hartwig.hmftools.common.exception.EmptyFolderException;
import com.hartwig.hmftools.common.exception.FolderDoesNotExistException;
import com.hartwig.hmftools.common.exception.HartwigException;

import org.jetbrains.annotations.NotNull;

@FunctionalInterface
public interface FolderChecker {

    @NotNull
    String checkFolder(@NotNull String directory) throws IOException, HartwigException;

    @NotNull
    static FolderChecker build() {
        return (directory) -> {
            final File folder = new File(directory);
            checkIfDirectoryExist(folder);
            checkIfIsDirectory(folder);
            checkIfDirectoryIsNotEmpty(folder);
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

    static void checkIfDirectoryIsNotEmpty(@NotNull final File directory) throws IOException, EmptyFolderException {
        final Path path = directory.toPath();
        if (isDirectoryEmpty(path)) {
            throw new EmptyFolderException(path.toString());
        }
    }

    static boolean isDirectoryEmpty(@NotNull final Path path) throws IOException {
        try (final DirectoryStream<Path> directoryStream = Files.newDirectoryStream(path, getHiddenFilesFilter())) {
            return !directoryStream.iterator().hasNext();
        }
    }

    @NotNull
    static DirectoryStream.Filter<Path> getHiddenFilesFilter() {
        return entry -> !Files.isHidden(entry);
    }
}
