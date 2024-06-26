package com.hartwig.hmftools.bamtools.metrics;

import static java.lang.String.format;

import static com.hartwig.hmftools.bamtools.common.CommonUtils.BT_LOGGER;
import static com.hartwig.hmftools.common.metrics.BamMetricsSummary.BAM_METRICS_FILE_ID;
import static com.hartwig.hmftools.common.utils.file.CommonFields.FLD_CHROMOSOME;
import static com.hartwig.hmftools.common.utils.file.CommonFields.FLD_POSITION_END;
import static com.hartwig.hmftools.common.utils.file.CommonFields.FLD_POSITION_START;
import static com.hartwig.hmftools.common.utils.file.FileDelimiters.TSV_DELIM;
import static com.hartwig.hmftools.common.utils.file.FileWriterUtils.createBufferedWriter;

import java.io.BufferedWriter;
import java.io.IOException;
import java.util.List;
import java.util.StringJoiner;

import com.hartwig.hmftools.common.region.ChrBaseRegion;

public class TargetRegionStats extends ChrBaseRegion
{
    public final ReadCounts FragmentCounts;

    public TargetRegionStats(final ChrBaseRegion region)
    {
        super(region.Chromosome, region.start(), region.end());
        FragmentCounts = new ReadCounts();
    }

    public String toString()
    {
        return format("%s fragments(%s)", super.toString(), FragmentCounts);
    }

    public static BufferedWriter initialiseWriter(final MetricsConfig config)
    {
        try
        {
            String filename = config.OutputDir + config.SampleId + BAM_METRICS_FILE_ID + ".target_regions.tsv";

            BufferedWriter writer = createBufferedWriter(filename, false);

            StringJoiner header = new StringJoiner(TSV_DELIM);

            header.add(FLD_CHROMOSOME);
            header.add(FLD_POSITION_START);
            header.add(FLD_POSITION_END);
            header.add("TotalFragments");
            header.add("DuplicateFragments");
            header.add("DualStrandFragments");

            writer.write(header.toString());
            writer.newLine();
            return writer;
        }
        catch(IOException e)
        {
            BT_LOGGER.error("failed to initialise target regions writer: {}", e.toString());
            return null;
        }
    }

    public synchronized static void writeStatistics(final BufferedWriter writer, final List<TargetRegionStats> targetRegionMetrics)
    {
        if(writer == null)
            return;

        try
        {
            for(TargetRegionStats regionStats : targetRegionMetrics)
            {
                StringJoiner sj = new StringJoiner(TSV_DELIM);
                sj.add(regionStats.Chromosome);
                sj.add(String.valueOf(regionStats.start()));
                sj.add(String.valueOf(regionStats.end()));
                sj.add(String.valueOf(regionStats.FragmentCounts.Total));
                sj.add(String.valueOf(regionStats.FragmentCounts.Duplicates));
                sj.add(String.valueOf(regionStats.FragmentCounts.DualStrand));

                writer.write(sj.toString());
                writer.newLine();
            }
        }
        catch(IOException e)
        {
            BT_LOGGER.error("failed to initialise target regions writer: {}", e.toString());
            return;
        }
    }
}
