package com.hartwig.hmftools.common.variant.vcf;

import java.io.BufferedReader;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.function.Predicate;

import com.hartwig.hmftools.common.exception.HartwigException;
import com.hartwig.hmftools.common.io.path.PathExtensionFinder;
import com.hartwig.hmftools.common.variant.GermlineVariant;
import com.hartwig.hmftools.common.variant.GermlineVariantFactory;

import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

public final class VCFFileStreamer {

    private static final Predicate<String> VCF_DATA_LINE_PREDICATE = new VCFDataLinePredicate();

    private VCFFileStreamer() {
    }

    @NotNull
    public static BufferedReader getVCFReader(@NotNull final String basePath, @NotNull final String fileExtension)
            throws IOException, HartwigException {
        final Path vcfPath = PathExtensionFinder.build().findPath(basePath, fileExtension);
        return Files.newBufferedReader(vcfPath);
    }

    @Nullable
    public static GermlineVariant nextVariant(@NotNull final BufferedReader reader) throws IOException {
        String line = safeReadLine(reader);
        while (line != null) {
            if (VCF_DATA_LINE_PREDICATE.test(line)) {
                final GermlineVariant potentialVariant = GermlineVariantFactory.fromVCFLine(line);
                if (potentialVariant != null) {
                    return potentialVariant;
                } else {
                    line = safeReadLine(reader);
                }
            } else {
                line = safeReadLine(reader);
            }
        }
        reader.close();
        return null;
    }

    @Nullable
    private static String safeReadLine(@NotNull final BufferedReader reader) {
        try {
            return reader.readLine();
        } catch (IOException e) {
            return null;
        }
    }
}
