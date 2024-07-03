package com.hartwig.hmftools.purple.fittingsnv;

import static java.lang.String.format;

import static com.hartwig.hmftools.common.utils.file.FileDelimiters.TSV_DELIM;
import static com.hartwig.hmftools.common.utils.file.FileWriterUtils.checkAddDirSeparator;

import java.io.File;
import java.io.IOException;
import java.nio.file.Files;
import java.util.List;
import java.util.StringJoiner;

import com.google.common.collect.Lists;

public final class PeakModelFile
{
    private static final String EXTENSION = ".purple.somatic.clonality.tsv";

    public static String generateFilename(final String basePath, final String sample)
    {
        return checkAddDirSeparator(basePath) + sample + EXTENSION;
    }

    public static void write(final String filename, final List<PeakModelData> model) throws IOException
    {
        Files.write(new File(filename).toPath(), toLines(model));
    }

    private static List<String> toLines(final List<PeakModelData> model)
    {
        final List<String> lines = Lists.newArrayList();
        lines.add(header());
        model.stream().map(PeakModelFile::toString).forEach(lines::add);
        return lines;
    }

    private static String header()
    {
        return new StringJoiner(TSV_DELIM, "", "")
                .add("peak")
                .add("bucket")
                .add("bucketWeight")
                .add("peakAvgWeight")
                .add("isValid")
                .add("isSubclonal")
                .toString();
    }

    private static String toString(final PeakModelData ratio)
    {
        return new StringJoiner(TSV_DELIM)
                .add(format("%.4f", ratio.Peak))
                .add(format("%.4f", ratio.Bucket))
                .add(format("%.4f", ratio.BucketWeight))
                .add(format("%.4f", ratio.PeakAvgWeight))
                .add(String.valueOf(ratio.IsValid))
                .add(String.valueOf(ratio.IsSubclonal))
                .toString();
    }
}
