package com.hartwig.hmftools.redux.consensus;

import static java.lang.String.format;

import static com.hartwig.hmftools.common.utils.file.FileDelimiters.TSV_DELIM;
import static com.hartwig.hmftools.common.utils.file.FileWriterUtils.createBufferedWriter;
import static com.hartwig.hmftools.redux.ReduxConfig.RD_LOGGER;

import java.io.BufferedWriter;
import java.io.IOException;
import java.util.Map;
import java.util.StringJoiner;

import com.google.common.collect.Maps;
import com.hartwig.hmftools.redux.ReduxConfig;

public class UltimaStats
{
    public final Map<Integer,Integer> LowQualCountFrequency;

    public UltimaStats()
    {
        LowQualCountFrequency = Maps.newHashMap();
    }

    public void addLowQualCount(int lowQualBaseCount)
    {
        addLowQualCount(lowQualBaseCount, 1);
    }

    public void addLowQualCount(int lowQualBaseCount, int count)
    {
        LowQualCountFrequency.put(lowQualBaseCount, LowQualCountFrequency.getOrDefault(lowQualBaseCount, 0) + count);
    }

    public void merge(final UltimaStats other)
    {
        for(Map.Entry<Integer,Integer> entry : other.LowQualCountFrequency.entrySet())
        {
            Integer count = LowQualCountFrequency.getOrDefault(entry.getKey(), 0);
            LowQualCountFrequency.put(entry.getKey(), entry.getValue() + count);
        }
    }

    public void writeStats(final ReduxConfig config)
    {
        try
        {
            String filename = config.formFilename("ultima_stats");
            BufferedWriter writer = createBufferedWriter(filename, false);

            StringJoiner sj = new StringJoiner(TSV_DELIM);
            writer.write("LowQualBases\tCount");
            writer.newLine();

            for(Map.Entry<Integer,Integer> entry : LowQualCountFrequency.entrySet())
            {
                writer.write(format("%d\t%d", entry.getKey(), entry.getValue()));
                writer.newLine();
            }

            writer.close();
        }
        catch(IOException e)
        {
            RD_LOGGER.error(" failed to write Ultima stats: {}", e.toString());
        }
    }
}
