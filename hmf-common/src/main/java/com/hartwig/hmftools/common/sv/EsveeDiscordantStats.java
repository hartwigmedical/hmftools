package com.hartwig.hmftools.common.sv;

import static com.hartwig.hmftools.common.sv.EsveeCommon.ESVEE_FILE_ID;
import static com.hartwig.hmftools.common.utils.file.FileDelimiters.TSV_DELIM;
import static com.hartwig.hmftools.common.utils.file.FileDelimiters.TSV_EXTENSION;
import static com.hartwig.hmftools.common.utils.file.FileWriterUtils.createBufferedWriter;

import java.io.BufferedWriter;
import java.io.File;
import java.io.IOException;
import java.nio.file.Files;
import java.util.List;
import java.util.StringJoiner;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;

public class EsveeDiscordantStats
{
    public final long TotalReads; // ie from the full original BAM(s)
    public final long PrepReads; // written by prep, so candidate discordant reads (not fragments)

    // note these counts are in fragment terms hence see doubling in rate calc
    public final long[] TypeCounts;

    public EsveeDiscordantStats(final long totalReads, final long prepReads, final long[] typeCounts)
    {
        TotalReads = totalReads;
        PrepReads = prepReads;
        TypeCounts = typeCounts;
    }

    public static final String PREP_DISC_STATS_FILE_ID = "disc_stats" + TSV_EXTENSION;

    private static String FILE_EXTENSION = ESVEE_FILE_ID + "." + PREP_DISC_STATS_FILE_ID;

    public static String generateFilename(String basePath, String sample)
    {
        return basePath + File.separator + sample + FILE_EXTENSION;
    }

    private static final String FLD_TOTAL_READS = "TotalReads";
    private static final String FLD_PREP_READS = "PrepReads";

    private static final Logger LOGGER = LogManager.getLogger(EsveeDiscordantStats.class);

    public static void write(final String filename, EsveeDiscordantStats stats)
    {
        try
        {
            BufferedWriter writer = createBufferedWriter(filename, false);

            StringJoiner sj = new StringJoiner(TSV_DELIM);
            sj.add(FLD_TOTAL_READS);
            sj.add(FLD_PREP_READS);

            for(DiscordantFragType type : DiscordantFragType.values())
            {
                sj.add(type.toString());
            }

            writer.write(sj.toString());
            writer.newLine();

            sj = new StringJoiner(TSV_DELIM);
            sj.add(String.valueOf(stats.TotalReads));
            sj.add(String.valueOf(stats.PrepReads));

            for(DiscordantFragType type : DiscordantFragType.values())
            {
                sj.add(String.valueOf(stats.TypeCounts[type.ordinal()]));
            }

            writer.write(sj.toString());
            writer.newLine();

            writer.close();
        }
        catch(IOException e)
        {
            LOGGER.error("failed to write Esvee discordant stats length file: {}", e.toString());
        }
    }

    public static EsveeDiscordantStats read(final String filename)
    {
        try
        {
            List<String> lines = Files.readAllLines(new File(filename).toPath());

            String[] values = lines.get(1).split(TSV_DELIM);

            long totalReads = Long.parseLong(values[0]);
            long prepReads = 0;

            long[] typeCounts = new long[DiscordantFragType.values().length];

            prepReads = Long.parseLong(values[1]);

            for(int i = 0; i < typeCounts.length; ++i)
            {
                typeCounts[i] = Long.parseLong(values[i + 2]);
            }

            return new EsveeDiscordantStats(totalReads, prepReads, typeCounts);
        }
        catch(Exception e)
        {
            LOGGER.error("failed to load Esvee discordant stat file({}): {}", filename, e.toString());
            return null;
        }
    }
}
