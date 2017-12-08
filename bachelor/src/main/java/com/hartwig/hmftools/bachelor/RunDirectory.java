package com.hartwig.hmftools.bachelor;

import java.io.File;
import java.io.IOException;
import java.nio.file.FileVisitOption;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.stream.Stream;

import org.jetbrains.annotations.Nullable;

class RunDirectory {

    final Path prefix;
    final File germline;
    final File somatic;
    final File copyNumber;
    final File structuralVariants;

    RunDirectory(final Path runDirectory) {
        prefix = runDirectory;
        germline = findGermline();
        somatic = findSomatic();
        copyNumber = findCopyNumber();
        structuralVariants = findStructuralVariants();
    }

    @Nullable
    private File findGermline() {
        try {
            try (final Stream<Path> stream = Files.walk(prefix.toRealPath(), 1, FileVisitOption.FOLLOW_LINKS)) {
                return stream.filter(p -> p.toString().endsWith("GoNLv5.vcf") || p.toString().endsWith("annotated.vcf"))
                        .map(Path::toFile)
                        .findFirst()
                        .orElse(null);
            }
        } catch (final IOException e) {
            return null;
        }
    }

    @Nullable
    private File findSomatic() {
        try {
            try (final Stream<Path> stream = Files.walk(prefix, FileVisitOption.FOLLOW_LINKS)) {
                return stream.filter(p -> p.toString().endsWith("_post_processed.vcf") || p.toString().endsWith("_melted.vcf"))
                        .map(Path::toFile)
                        .findFirst()
                        .orElse(null);
            }
        } catch (final IOException e) {
            return null;
        }
    }

    @Nullable
    private File findCopyNumber() {
        try {
            try (final Stream<Path> stream = Files.walk(prefix, FileVisitOption.FOLLOW_LINKS)) {
                return stream.filter(p -> p.toString().endsWith("purple.cnv")).map(Path::toFile).findFirst().orElse(null);
            }
        } catch (final IOException e) {
            return null;
        }
    }

    @Nullable
    private File findStructuralVariants() {
        try {
            try (final Stream<Path> stream = Files.walk(prefix, FileVisitOption.FOLLOW_LINKS)) {
                return stream.filter(p -> p.toString().endsWith("bpi.vcf")).map(Path::toFile).findFirst().orElse(null);
            }
        } catch (final IOException e) {
            return null;
        }
    }

    String getPatientID() {
        final String[] split = prefix.getFileName().toString().split("_");
        return split[split.length - 1];
    }
}
