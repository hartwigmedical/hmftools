package com.hartwig.hmftools.serve;

import org.jetbrains.annotations.NotNull;

public enum RefGenomeVersion {
    V37("37"),
    V38("38");

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
        if (!filePath.contains(".")) {
            throw new IllegalStateException("Cannot include ref genome version in file path that has no extension: " + filePath);
        }

        int extensionStart = filePath.lastIndexOf(".");
        return filePath.substring(0, extensionStart) + "." + this.identifier + filePath.substring(extensionStart);
    }
}
