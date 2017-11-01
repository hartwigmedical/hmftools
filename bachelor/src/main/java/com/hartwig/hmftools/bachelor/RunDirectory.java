package com.hartwig.hmftools.bachelor;

import java.io.File;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Path;

import org.jetbrains.annotations.Nullable;

class RunDirectory {

    final Path prefix;

    RunDirectory(final Path runDirectory) {
        prefix = runDirectory;
    }

    @Nullable
    File findGermline() {
        try {
            return Files.walk(prefix, 1)
                    .filter(p -> p.toString().endsWith("GoNLv5.vcf") || p.toString().endsWith("annotated.vcf"))
                    .map(Path::toFile)
                    .findFirst()
                    .orElse(null);
        } catch (final IOException e) {
            return null;
        }
    }

    @Nullable
    File findSomatic() {
        try {
            return Files.walk(prefix)
                    .filter(p -> p.toString().endsWith("_post_processed.vcf") || p.toString().endsWith("_melted.vcf"))
                    .map(Path::toFile)
                    .findFirst()
                    .orElse(null);
        } catch (final IOException e) {
            return null;
        }
    }

    @Nullable
    File findCopyNumber() {
        try {
            return Files.walk(prefix).filter(p -> p.toString().endsWith("purple.cnv")).map(Path::toFile).findFirst().orElse(null);
        } catch (final IOException e) {
            return null;
        }
    }

    @Nullable
    File findStructuralVariants() {
        try {
            return Files.walk(prefix).filter(p -> p.toString().endsWith("bpi.vcf")).map(Path::toFile).findFirst().orElse(null);
        } catch (final IOException e) {
            return null;
        }
    }
}
