package com.hartwig.hmftools.common.serve;

import org.jetbrains.annotations.NotNull;

public enum RefGenomeVersion {
    V37("37"),
    V38("38");

    private static final String GZIP_EXTENSION = ".gz";

    @NotNull
    private final String identifier;

    RefGenomeVersion(@NotNull final String identifier) {
        this.identifier = identifier;
    }

    @NotNull
    public static RefGenomeVersion fromIdentifier(@NotNull String identifier) {
        for (RefGenomeVersion version : values()) {
            if (version.identifier.equals(identifier)) {
                return version;
            }
        }

        throw new IllegalStateException("Could not find ref genome version with identifier '" + identifier + "'");
    }

    @NotNull
    public String addVersionToFilePath(@NotNull String filePath) {
        String modifiedFilePath = filePath;
        if (filePath.endsWith(GZIP_EXTENSION)) {
            modifiedFilePath = filePath.substring(0, filePath.indexOf(GZIP_EXTENSION));
        }

        if (!modifiedFilePath.contains(".")) {
            throw new IllegalStateException("Cannot include ref genome version in file path that has no proper extension: " + filePath);
        }

        int extensionStart = modifiedFilePath.lastIndexOf(".");
        String versionedFilePath =
                modifiedFilePath.substring(0, extensionStart) + "." + this.identifier + modifiedFilePath.substring(extensionStart);

        if (filePath.endsWith(GZIP_EXTENSION)) {
            versionedFilePath = versionedFilePath + GZIP_EXTENSION;
        }

        return versionedFilePath;
    }
}
