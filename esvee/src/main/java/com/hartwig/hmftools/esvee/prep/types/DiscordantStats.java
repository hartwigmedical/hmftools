package com.hartwig.hmftools.esvee.prep.types;

import static java.lang.Math.abs;
import static java.lang.String.format;

import static com.hartwig.hmftools.common.utils.file.FileDelimiters.ITEM_DELIM;
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
import java.lang.reflect.Type;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.Arrays;
import java.util.List;
import java.util.StringJoiner;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.esvee.common.FragmentLengthBounds;
import com.hartwig.hmftools.esvee.prep.PrepConfig;

public class DiscordantStats
{
    public final long ProcessedReads;
    public final long WrittenReads;

    // note these counts are in fragment terms hence see doubling in rate calc
    public final long[] TypeCounts;

    public static final int SHORT_INV_LENGTH = 5000;
    public static final int SHORT_LENGTH = 1000;
    public static final int MEDIUM_LENGTH = 10000;
    public static final int LONG_LENGTH = 100000;

    private enum DiscordantType
    {
        Translocation,
        ShortInv,
        DelDup10K,
        DelDup100K,
        DelDupLong,
        Inv100K,
        InvLong;
    }

    public DiscordantStats(final long processedReads, final long writtenReads, final long[] typeCounts)
    {
        ProcessedReads = processedReads;
        WrittenReads = writtenReads;
        TypeCounts = typeCounts;
    }

    public DiscordantStats()
    {
        ProcessedReads = 0;
        WrittenReads = 0;
        TypeCounts = new long[DiscordantType.values().length];
    }

    public void addRead(final PrepRead read)
    {
        if(!read.Chromosome.equals(read.MateChromosome))
        {
            ++TypeCounts[DiscordantType.Translocation.ordinal()];
            return;
        }

        int distance = abs(read.record().getInferredInsertSize());

        if(read.orientation() == read.mateOrientation())
        {
            if(distance < SHORT_INV_LENGTH)
                TypeCounts[DiscordantType.ShortInv.ordinal()] += 2; // assumes mate is in the same read group so not also registered
            else if(distance < LONG_LENGTH)
                ++TypeCounts[DiscordantType.Inv100K.ordinal()];
            else
                ++TypeCounts[DiscordantType.InvLong.ordinal()];
        }
        else
        {
            if(distance < MEDIUM_LENGTH)
                ++TypeCounts[DiscordantType.DelDup10K.ordinal()];
            else if(distance < LONG_LENGTH)
                ++TypeCounts[DiscordantType.DelDup100K.ordinal()];
            else
                ++TypeCounts[DiscordantType.DelDupLong.ordinal()];
        }
    }

    public void add(final DiscordantStats stats)
    {
        for(int i = 0; i < TypeCounts.length; ++i)
        {
            TypeCounts[i] += stats.TypeCounts[i];
        }
    }

    public double shortInversionRate() { return ProcessedReads > 0 ? TypeCounts[DiscordantType.ShortInv.ordinal()] / (double)ProcessedReads : 0; }

    public double discordantRate()
    {
        double totalDiscordant = Arrays.stream(TypeCounts).sum();
        return ProcessedReads > 0 ? totalDiscordant / (double)ProcessedReads : 0;
    }

    public String toString()
    {
        StringJoiner sj = new StringJoiner(ITEM_DELIM);

        for(DiscordantType type : DiscordantType.values())
        {
            sj.add(format("%s=%d", type, TypeCounts[type.ordinal()]));
        }

        return format("totalReads=%d writtenReads=%d %s", ProcessedReads, WrittenReads, sj);
    }

    private static final String FLD_TOTAL_READS = "TotalReads";
    private static final String FLD_WRITTEN_READS = "WrittenReads";

    public static void writeDiscordantStats(
            final PrepConfig config, final long processedReads, final long writtenReads, final DiscordantStats discordantStats)
    {
        try
        {
            final String outputFileName = config.formFilename(DISCORDANT_STATS);

            BufferedWriter writer = createBufferedWriter(outputFileName, false);

            StringJoiner sj = new StringJoiner(TSV_DELIM);
            sj.add(FLD_TOTAL_READS);
            sj.add(FLD_WRITTEN_READS);

            for(DiscordantType type : DiscordantType.values())
            {
                sj.add(type.toString());
            }

            writer.write(sj.toString());
            writer.newLine();

            sj = new StringJoiner(TSV_DELIM);
            sj.add(String.valueOf(processedReads));
            sj.add(String.valueOf(writtenReads));

            for(DiscordantType type : DiscordantType.values())
            {
                sj.add(String.valueOf(discordantStats.TypeCounts[type.ordinal()]));
            }

            writer.write(sj.toString());
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

            long processedReads = Long.parseLong(values[0]);
            long writtenReads = 0;

            long[] typeCounts = new long[DiscordantType.values().length];

            // backwards compatibility
            if(values.length == 4)
            {
                typeCounts[DiscordantType.ShortInv.ordinal()] = Long.parseLong(values[1]);
                typeCounts[DiscordantType.Translocation.ordinal()] = Long.parseLong(values[2]);
                typeCounts[DiscordantType.DelDupLong.ordinal()] = Long.parseLong(values[3]);
            }
            else
            {
                writtenReads = Long.parseLong(values[1]);

                for(int i = 0; i < typeCounts.length; ++i)
                {
                    typeCounts[i] = Long.parseLong(values[i + 2]);
                }
            }

            return new DiscordantStats(processedReads, writtenReads, typeCounts);
        }
        catch(Exception e)
        {
            SV_LOGGER.warn("failed to read discordant read statistics file: {}", e.toString());
            return new DiscordantStats();
        }
    }
}
