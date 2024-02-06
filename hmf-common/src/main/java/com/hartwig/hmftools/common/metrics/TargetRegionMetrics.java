package com.hartwig.hmftools.common.metrics;

import static java.lang.String.format;

import static com.hartwig.hmftools.common.metrics.BamMetricsSummary.BAM_METRICS_FILE_ID;
import static com.hartwig.hmftools.common.metrics.BamMetricsSummary.DUAL_STRAND_COLUMN;
import static com.hartwig.hmftools.common.metrics.BamMetricsSummary.DUPLICATES_COLUMN;
import static com.hartwig.hmftools.common.metrics.BamMetricsSummary.TOTAL_READS_COLUMN;
import static com.hartwig.hmftools.common.utils.file.CommonFields.FLD_CHROMOSOME;
import static com.hartwig.hmftools.common.utils.file.CommonFields.FLD_POSITION_END;
import static com.hartwig.hmftools.common.utils.file.CommonFields.FLD_POSITION_START;
import static com.hartwig.hmftools.common.utils.file.FileDelimiters.TSV_DELIM;
import static com.hartwig.hmftools.common.utils.file.FileReaderUtils.createFieldsIndexMap;
import static com.hartwig.hmftools.common.utils.file.FileWriterUtils.checkAddDirSeparator;

import java.io.File;
import java.io.IOException;
import java.nio.file.Files;
import java.util.List;
import java.util.Map;
import java.util.StringJoiner;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.region.ChrBaseRegion;

import org.immutables.value.Value;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

@Value.Immutable
@Value.Style(passAnnotations = { NotNull.class, Nullable.class })
public abstract class TargetRegionMetrics
{
    public abstract ChrBaseRegion region();
    public abstract long totalReads();
    public abstract long duplicateReads();
    public abstract long dualStrandReads();

    private static final String FILE_EXTENSION = ".target_regions.tsv";

    public static String generateFilename(final String basePath, final String sampleId)
    {
        return checkAddDirSeparator(basePath) + sampleId + BAM_METRICS_FILE_ID + FILE_EXTENSION;
    }

    public static void write(final String filename, final List<TargetRegionMetrics> targetRegions) throws IOException
    {
        List<String> lines = Lists.newArrayList();

        StringJoiner header = new StringJoiner(TSV_DELIM);

        header.add(FLD_CHROMOSOME);
        header.add(FLD_POSITION_START);
        header.add(FLD_POSITION_END);
        header.add(TOTAL_READS_COLUMN);
        header.add(DUPLICATES_COLUMN);
        header.add(DUAL_STRAND_COLUMN);

        lines.add(header.toString());

        for(TargetRegionMetrics region : targetRegions)
        {
            StringJoiner values = new StringJoiner(TSV_DELIM);
            values.add(region.region().Chromosome);
            values.add(String.valueOf(region.region().start()));
            values.add(String.valueOf(region.region().end()));
            values.add(String.valueOf(region.totalReads()));
            values.add(String.valueOf(region.duplicateReads()));
            values.add(String.valueOf(region.dualStrandReads()));

            lines.add(values.toString());
        }

        Files.write(new File(filename).toPath(), lines);
    }

    public static List<TargetRegionMetrics> read(final String filename) throws IOException
    {
        List<String> lines = Files.readAllLines(new File(filename).toPath());

        String header = lines.get(0);
        lines.remove(0);

        Map<String,Integer> fieldsIndexMap = createFieldsIndexMap(header, TSV_DELIM);

        List<TargetRegionMetrics> targetRegions = Lists.newArrayList();

        for(String line : lines)
        {
            String[] values = line.split(TSV_DELIM, -1);

            ChrBaseRegion region = new ChrBaseRegion(
                    values[fieldsIndexMap.get(FLD_CHROMOSOME)],
                    Integer.parseInt(values[fieldsIndexMap.get(FLD_POSITION_START)]),
                    Integer.parseInt(values[fieldsIndexMap.get(FLD_POSITION_END)]));

            targetRegions.add(ImmutableTargetRegionMetrics.builder()
                    .region(region)
                    .totalReads(Integer.parseInt(values[fieldsIndexMap.get(TOTAL_READS_COLUMN)]))
                    .duplicateReads(Integer.parseInt(values[fieldsIndexMap.get(DUPLICATES_COLUMN)]))
                    .dualStrandReads(Integer.parseInt(values[fieldsIndexMap.get(DUAL_STRAND_COLUMN)]))
                    .build());
        }

        return targetRegions;
    }
}
