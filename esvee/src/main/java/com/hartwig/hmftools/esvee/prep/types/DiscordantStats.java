package com.hartwig.hmftools.esvee.prep.types;

import static java.lang.String.format;

import static com.hartwig.hmftools.common.utils.file.FileWriterUtils.createBufferedWriter;
import static com.hartwig.hmftools.esvee.assembly.AssemblyConfig.SV_LOGGER;
import static com.hartwig.hmftools.esvee.prep.types.WriteType.DISCORDANT_STATS;
import static com.hartwig.hmftools.esvee.prep.types.WriteType.FRAGMENT_LENGTH_DIST;

import java.io.BufferedWriter;
import java.io.IOException;
import java.util.StringJoiner;

import com.hartwig.hmftools.esvee.prep.PrepConfig;

public class DiscordantStats
{
    public int ShortInversion;
    public int Translocation;
    public int Other;

    public static final int SHORT_INV_LENGTH = 5000;

    public DiscordantStats()
    {
        ShortInversion = 0;
        Translocation = 0;
        Other = 0;
    }

    public void add(final DiscordantStats stats)
    {
        ShortInversion += stats.ShortInversion;
        Translocation += stats.Translocation;
        Other += stats.Other;
    }

    public String toString() { return format("shortInv(%d) bnd(%d) other(%d)", ShortInversion, Translocation, Other); }

    public static void writeDiscordantStats(final PrepConfig config, final long totalReads, final DiscordantStats discordantStats)
    {
        try
        {
            final String outputFileName = config.formFilename(DISCORDANT_STATS);

            BufferedWriter writer = createBufferedWriter(outputFileName, false);

            writer.write("TotalReads\tShortInversions\tTranslocations\tOther");
            writer.newLine();
            writer.write(String.format("%d\t%d\t%d\t%d",
                    totalReads, discordantStats.ShortInversion, discordantStats.Translocation, discordantStats.Other));
            writer.newLine();

            writer.close();
        }
        catch(IOException e)
        {
            SV_LOGGER.error("failed to write fragment length file: {}", e.toString());
        }
    }
}
