package com.hartwig.hmftools.common.variant.vcfloader;

import java.io.IOException;
import java.nio.file.Path;
import java.util.List;
import java.util.stream.Collectors;

import com.hartwig.hmftools.common.exception.HealthChecksException;
import com.hartwig.hmftools.common.io.path.PathExtensionFinder;
import com.hartwig.hmftools.common.io.reader.LineReader;
import com.hartwig.hmftools.common.variant.SomaticVariant;
import com.hartwig.hmftools.common.variant.SomaticVariantFactory;
import com.hartwig.hmftools.common.variant.predicate.VCFPassDataLinePredicate;

import org.jetbrains.annotations.NotNull;

public final class VCFFileLoader {

    private VCFFileLoader() {
    }

    @NotNull
    public static List<SomaticVariant> loadSomaticVCF(@NotNull String runDirectory, @NotNull String fileExtension)
            throws IOException, HealthChecksException {
        final Path vcfPath = PathExtensionFinder.build().findPath(runDirectory, fileExtension);
        final List<String> lines = LineReader.build().readLines(vcfPath, new VCFPassDataLinePredicate());
        return lines.stream().map(SomaticVariantFactory::fromVCFLine).collect(Collectors.toList());
    }
}
