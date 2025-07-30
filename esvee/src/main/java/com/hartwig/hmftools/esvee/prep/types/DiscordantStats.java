package com.hartwig.hmftools.esvee.prep.types;

import static java.lang.Math.abs;
import static java.lang.String.format;

import static com.hartwig.hmftools.common.bam.SamRecordUtils.inferredInsertSizeAbs;
import static com.hartwig.hmftools.common.utils.file.FileDelimiters.ITEM_DELIM;
import static com.hartwig.hmftools.common.utils.file.FileDelimiters.TSV_DELIM;
import static com.hartwig.hmftools.common.utils.file.FileWriterUtils.createBufferedWriter;
import static com.hartwig.hmftools.esvee.assembly.AssemblyConfig.SV_LOGGER;
import static com.hartwig.hmftools.esvee.prep.types.WriteType.DISCORDANT_STATS;

import java.io.BufferedWriter;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.Arrays;
import java.util.List;
import java.util.StringJoiner;

import com.hartwig.hmftools.common.genome.region.Orientation;
import com.hartwig.hmftools.esvee.prep.PrepConfig;

public class DiscordantStats
{
    public final long TotalReads; // ie from the full original BAM(s)
    public final long PrepReads; // written by prep, so candidate discordant reads (not fragments)

    // note these counts are in fragment terms hence see doubling in rate calc
    public final long[] TypeCounts;

    public DiscordantStats(final long totalReads, final long prepReads, final long[] typeCounts)
    {
        TotalReads = totalReads;
        PrepReads = prepReads;
        TypeCounts = typeCounts;
    }

    public DiscordantStats()
    {
        TotalReads = 0;
        PrepReads = 0;
        TypeCounts = new long[DiscordantFragType.values().length];
    }

    private static final int LENGTH_SHORT = 1000;
    private static final int LENGTH_5K = 5000;
    private static final int LENGTH_100K = 100_000;

    public void processReadGroup(final ReadGroup readGroup)
    {
        PrepRead read = readGroup.reads().stream().filter(x -> !x.isSupplementaryAlignment()).findFirst().orElse(null);

        addToCounts(
                read.Chromosome, read.MateChromosome, read.AlignmentStart, read.MatePosStart, read.orientation(), read.mateOrientation(),
                inferredInsertSizeAbs(read.record()));
    }

    private void addToCounts(
            final String chromosome1, final String chromosome2, int posStart1, int posStart2,
            final Orientation orientation1, final Orientation orientation2, int insertSize)
    {
        if(!chromosome1.equals(chromosome2))
        {
            ++TypeCounts[DiscordantFragType.Translocation.ordinal()];
            return;
        }

        int distance = insertSize;

        if(orientation1 == orientation2)
        {
            if(distance < LENGTH_SHORT)
                TypeCounts[DiscordantFragType.InvLt1K.ordinal()] += 2; // assumes mate is in the same read group so not also registered
            else if(distance < LENGTH_5K)
                ++TypeCounts[DiscordantFragType.Inv1To5K.ordinal()];
            else if(distance < LENGTH_100K)
                ++TypeCounts[DiscordantFragType.Inv5To100K.ordinal()];
            else
                ++TypeCounts[DiscordantFragType.InvGt100K.ordinal()];
        }
        else
        {
            if((posStart1 < posStart2) == orientation1.isForward())
            {
                // DEL-type orientation
                if(distance < LENGTH_5K)
                    ++TypeCounts[DiscordantFragType.Del1To5K.ordinal()];
                else if(distance < LENGTH_100K)
                    ++TypeCounts[DiscordantFragType.Del5To100K.ordinal()];
                else
                    ++TypeCounts[DiscordantFragType.DelGt100K.ordinal()];
            }
            else
            {
                if(distance < LENGTH_5K)
                    ++TypeCounts[DiscordantFragType.Dup1To5K.ordinal()];
                else if(distance < LENGTH_100K)
                    ++TypeCounts[DiscordantFragType.Dup5To100K.ordinal()];
                else
                    ++TypeCounts[DiscordantFragType.DupGt100K.ordinal()];
            }
        }
    }

    public void add(final DiscordantStats stats)
    {
        for(int i = 0; i < TypeCounts.length; ++i)
        {
            TypeCounts[i] += stats.TypeCounts[i];
        }
    }

    public double shortInversionRate()
    {
        long shortInvCount = TypeCounts[DiscordantFragType.InvLt1K.ordinal()] + TypeCounts[DiscordantFragType.Inv1To5K.ordinal()];
        return TotalReads > 0 ? shortInvCount / (double) TotalReads : 0;
    }

    public double shortInversionFragmentRate()
    {
        long shortInvCount = TypeCounts[DiscordantFragType.InvLt1K.ordinal()];
        return TotalReads > 0 ? shortInvCount / (double) TotalReads : 0;
    }

    public double discordantRate()
    {
        double totalDiscordant = Arrays.stream(TypeCounts).sum();
        return TotalReads > 0 ? totalDiscordant / (double) TotalReads : 0;
    }

    public String toString()
    {
        StringJoiner sj = new StringJoiner(ITEM_DELIM);

        for(DiscordantFragType type : DiscordantFragType.values())
        {
            sj.add(format("%s=%d", type, TypeCounts[type.ordinal()]));
        }

        return format("totalReads=%d writtenReads=%d %s", TotalReads, PrepReads, sj);
    }

    private static final String FLD_TOTAL_READS = "TotalReads";
    private static final String FLD_PREP_READS = "PrepReads";

    public static void writeDiscordantStats(
            final PrepConfig config, final long totalReads, final long writtenReads, final DiscordantStats discordantStats)
    {
        try
        {
            final String outputFileName = config.formFilename(DISCORDANT_STATS);

            BufferedWriter writer = createBufferedWriter(outputFileName, false);

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
            sj.add(String.valueOf(totalReads));
            sj.add(String.valueOf(writtenReads));

            for(DiscordantFragType type : DiscordantFragType.values())
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

    public static DiscordantStats loadDiscordantStats(final String filename)
    {
        try
        {
            List<String> lines = Files.readAllLines(Paths.get(filename));
            String[] values = lines.get(1).split(TSV_DELIM);

            long totalReads = Long.parseLong(values[0]);
            long prepReads = 0;

            long[] typeCounts = new long[DiscordantFragType.values().length];

            prepReads = Long.parseLong(values[1]);

            for(int i = 0; i < typeCounts.length; ++i)
            {
                typeCounts[i] = Long.parseLong(values[i + 2]);
            }

            return new DiscordantStats(totalReads, prepReads, typeCounts);
        }
        catch(Exception e)
        {
            SV_LOGGER.warn("failed to read discordant read statistics file({}): {}", filename, e.toString());
            return new DiscordantStats();
        }
    }
}
