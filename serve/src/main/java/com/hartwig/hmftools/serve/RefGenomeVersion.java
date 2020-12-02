package com.hartwig.hmftools.serve;

import org.jetbrains.annotations.NotNull;

public enum RefGenomeVersion {
    HG19;

    @NotNull
    public String makeVersioned(@NotNull String filePath) {
        if (!filePath.contains(".")) {
            throw new IllegalStateException("Cannot include ref genome version in file path that has no extension: " + filePath);
        }

        int extensionStart = filePath.lastIndexOf(".");
        return filePath.substring(0, extensionStart) + "." + this.toString().toLowerCase() + filePath.substring(extensionStart);
    }
}
