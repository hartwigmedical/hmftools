package com.hartwig.hmftools.esvee.prep.types;

import static java.lang.String.format;

import static com.hartwig.hmftools.common.utils.file.FileDelimiters.TSV_DELIM;
import static com.hartwig.hmftools.common.utils.file.FileWriterUtils.createBufferedWriter;
import static com.hartwig.hmftools.esvee.assembly.AssemblyConfig.SV_LOGGER;
import static com.hartwig.hmftools.esvee.common.FileCommon.formPrepInputFilename;
import static com.hartwig.hmftools.esvee.common.FragmentLengthBounds.INVALID;
import static com.hartwig.hmftools.esvee.prep.PrepConstants.PREP_DISC_STATS_FILE_ID;
import static com.hartwig.hmftools.esvee.prep.PrepConstants.PREP_FRAG_LENGTH_FILE_ID;
import static com.hartwig.hmftools.esvee.prep.types.WriteType.DISCORDANT_STATS;
import static com.hartwig.hmftools.esvee.prep.types.WriteType.FRAGMENT_LENGTH_DIST;

import java.io.BufferedWriter;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.List;
import java.util.StringJoiner;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.esvee.common.FragmentLengthBounds;
import com.hartwig.hmftools.esvee.prep.PrepConfig;

public class DiscordantStats
{
    public final long TotalReads;

    // note these counts are in fragment terms hence see doubling in rate calc
    public long ShortInversion;
    public long Translocation;
    public long Other;

    public static final int SHORT_INV_LENGTH = 5000;

    public DiscordantStats(final long totalReads, final long shortInversion, final long translocation, final long other)
    {
        TotalReads = totalReads;
        ShortInversion = shortInversion;
        Translocation = translocation;
        Other = other;
    }

    public DiscordantStats()
    {
        TotalReads = 0;
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

    public double shortInversionRate() { return TotalReads > 0 ? 2.0 * ShortInversion / TotalReads : 0; }

    public String toString() { return format("totalReads(%d) shortInv(%d) bnd(%d) other(%d)", TotalReads, ShortInversion, Translocation, Other); }

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

    public static String formDiscordantStatsFilename(final String outputDir, final String sampleId)
    {
        return formPrepInputFilename(outputDir, sampleId, PREP_DISC_STATS_FILE_ID, null);
    }

    public static DiscordantStats loadDiscordantStats(final String filename)
    {
        try
        {
            List<String> lines = Files.readAllLines(Paths.get(filename));
            String[] values = lines.get(1).split(TSV_DELIM);

            return new DiscordantStats(
                    Long.parseLong(values[0]), Long.parseLong(values[1]), Long.parseLong(values[2]), Long.parseLong(values[3]));
        }
        catch(Exception e)
        {
            SV_LOGGER.warn("failed to read discordant read statsitics file: {}", e.toString());
            return new DiscordantStats();
        }
    }
}
