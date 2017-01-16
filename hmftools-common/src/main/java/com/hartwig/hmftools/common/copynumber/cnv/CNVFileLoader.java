package com.hartwig.hmftools.common.copynumber.cnv;

import java.io.IOException;
import java.nio.file.Path;
import java.util.List;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.copynumber.CopyNumber;
import com.hartwig.hmftools.common.copynumber.CopyNumberFactory;
import com.hartwig.hmftools.common.exception.HartwigException;
import com.hartwig.hmftools.common.io.path.PathPrefixSuffixFinder;
import com.hartwig.hmftools.common.io.reader.FileReader;

import org.jetbrains.annotations.NotNull;

public final class CNVFileLoader {

    private static final String COPYNUMBER_SUFFIX = ".bam_CNVs";

    private CNVFileLoader() {
    }

    @NotNull
    public static List<CopyNumber> loadCNV(@NotNull final String basePath, @NotNull final String sample,
            @NotNull final String fileExtension) throws IOException, HartwigException {
        final Path copynumberPath = PathPrefixSuffixFinder.build().findPath(basePath, sample, COPYNUMBER_SUFFIX);
        final List<String> lines = FileReader.build().readLines(copynumberPath);

        final List<CopyNumber> copyNumbers = Lists.newArrayList();
        for (final String line : lines) {
            copyNumbers.add(CopyNumberFactory.fromCNVLine(line));
        }
        return copyNumbers;
    }
}
