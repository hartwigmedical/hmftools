package com.hartwig.hmftools.common.variant.vcf;

import java.io.File;
import java.io.IOException;
import java.nio.file.Files;
import java.util.Collection;
import java.util.List;
import java.util.stream.Collectors;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.variant.SomaticVariant;
import com.hartwig.hmftools.common.variant.SomaticVariantFactory;

import org.jetbrains.annotations.NotNull;

public final class VCFFileWriter {

    private VCFFileWriter() {
    }

    public static void writeSomaticVCF(@NotNull final String filePath, @NotNull final VCFSomaticFile originalFile,
            @NotNull final List<SomaticVariant> newVariants) throws IOException {
        final Collection<String> vcfLines = Lists.newArrayList();
        vcfLines.addAll(originalFile.originalMetaInformationLines());
        vcfLines.add(originalFile.originalHeaderLine());
        vcfLines.addAll(newVariants.stream().map(SomaticVariantFactory::toVCFLine).collect(Collectors.toList()));

        Files.write(new File(filePath).toPath(), vcfLines);
    }
}
