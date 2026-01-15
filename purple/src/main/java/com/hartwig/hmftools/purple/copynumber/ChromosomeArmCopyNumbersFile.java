package com.hartwig.hmftools.purple.copynumber;

import static com.hartwig.hmftools.common.utils.file.FileWriterUtils.checkAddDirSeparator;

import java.io.File;
import java.io.IOException;
import java.nio.file.Files;
import java.util.Collection;
import java.util.List;

import com.google.common.collect.Lists;

public class ChromosomeArmCopyNumbersFile
{
    public static final String EXTENSION = ".purple.chromosome_arm.tsv";

    public static String generateFilename(final String basePath, final String sample)
    {
        return checkAddDirSeparator(basePath) + sample + EXTENSION;
    }

    public static void write(final String filename, final List<ChromosomeArmCopyNumber> segments) throws IOException
    {
        final Collection<String> lines = Lists.newArrayList();
        lines.add(ChromosomeArmCopyNumber.tsvFileHeader());
        segments.stream().map(ChromosomeArmCopyNumber::toTSV).forEach(lines::add);
        Files.write(new File(filename).toPath(), lines);
    }

    public static List<ChromosomeArmCopyNumber> read(final String filename) throws IOException
    {
        final List<String> lines = Files.readAllLines(new File(filename).toPath());
        final List<ChromosomeArmCopyNumber> result = Lists.newArrayList();
        for (int i = 1; i < lines.size(); i++) {
            result.add(ChromosomeArmCopyNumber.fromTsv(lines.get(i)));
        }
        return result;
    }
}
