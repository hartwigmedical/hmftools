package com.hartwig.hmftools.purple.targeted;

import static com.hartwig.hmftools.common.utils.file.FileWriterUtils.checkAddDirSeparator;

import java.io.File;
import java.io.IOException;
import java.nio.file.Files;
import java.util.Collection;
import java.util.List;

import com.google.common.collect.Lists;

public class TargetRegionsCopyNumberFile
{
    public static final String EXTENSION = ".purple.copynumbertargeted.tsv";

    public static String generateFilename(final String basePath, final String sample)
    {
        return checkAddDirSeparator(basePath) + sample + EXTENSION;
    }

    public static void write(final String filename, final List<TargetRegionsCopyNumber> segments) throws IOException
    {
        final Collection<String> lines = Lists.newArrayList();
        lines.add(TargetRegionsCopyNumber.tsvFileHeader());
        segments.stream().map(TargetRegionsCopyNumber::toTSV).forEach(lines::add);
        Files.write(new File(filename).toPath(), lines);
    }
}
