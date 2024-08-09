package com.hartwig.hmftools.orange.algo.purple;

import static com.hartwig.hmftools.common.utils.file.FileDelimiters.TSV_DELIM;
import static com.hartwig.hmftools.common.utils.file.FileReaderUtils.createFieldsIndexMap;
import static com.hartwig.hmftools.common.utils.file.FileWriterUtils.checkAddDirSeparator;

import java.io.File;
import java.io.IOException;
import java.nio.file.Files;
import java.util.List;
import java.util.Map;
import java.util.stream.Collectors;

import com.hartwig.hmftools.common.purple.GermlineStatus;

import org.jetbrains.annotations.NotNull;

public final class SegmentFile
{
    private static final String EXTENSION = ".purple.segment.tsv";

    @NotNull
    public static String generateFilename(final String basePath, final String sample)
    {
        return checkAddDirSeparator(basePath) + sample + EXTENSION;
    }

    @NotNull
    public static List<Segment> read(final String filePath) throws IOException
    {
        return fromLines(Files.readAllLines(new File(filePath).toPath()));
    }

    @NotNull
    static List<Segment> fromLines(final List<String> lines)
    {
        final Map<String, Integer> fieldsIndexMap = createFieldsIndexMap(lines.get(0), TSV_DELIM);
        lines.remove(0);

        int germlineStatusIndex = fieldsIndexMap.get("germlineStatus");
        int bafCountIndex = fieldsIndexMap.get("bafCount");
        int observedTumorRatioIndex = fieldsIndexMap.get("observedTumorRatio");

        return lines.stream().map(line ->
        {
            String[] values = line.split(TSV_DELIM, -1);
            return ImmutableSegment.builder()
                    .bafCount(Integer.parseInt(values[bafCountIndex]))
                    .observedTumorRatio(Double.parseDouble(values[observedTumorRatioIndex]))
                    .germlineStatus(GermlineStatus.valueOf(values[germlineStatusIndex]))
                    .build();
        }).collect(Collectors.toList());
    }
}