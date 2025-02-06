package com.hartwig.hmftools.esvee.prep.types;

import static java.lang.Math.abs;
import static java.lang.String.format;

import static com.hartwig.hmftools.common.bam.SamRecordUtils.inferredInsertSizeAbs;
import static com.hartwig.hmftools.common.region.BaseRegion.positionWithin;
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

import com.hartwig.hmftools.common.bam.SamRecordUtils;
import com.hartwig.hmftools.common.bam.SupplementaryReadData;
import com.hartwig.hmftools.common.genome.region.Orientation;
import com.hartwig.hmftools.esvee.prep.PrepConfig;

public class DiscordantStats
{
    public final long ProcessedReads;
    public final long PrepReads;

    // note these counts are in fragment terms hence see doubling in rate calc
    public final long[] TypeCounts;

    private enum DiscordantType
    {
        Translocation,
        InvLt1K,
        Inv1To5K,
        Inv5To100K,
        InvGt100K,
        Del1Lt1K,
        Del1To5K,
        Del5To100K,
        DelGt100K,
        Dup1To5K,
        Dup5To100K,
        DupGt100K;
    }

    public DiscordantStats(final long processedReads, final long prepReads, final long[] typeCounts)
    {
        ProcessedReads = processedReads;
        PrepReads = prepReads;
        TypeCounts = typeCounts;
    }

    public DiscordantStats()
    {
        ProcessedReads = 0;
        PrepReads = 0;
        TypeCounts = new long[DiscordantType.values().length];
    }

    private static final int LENGTH_SHORT = 1000;
    private static final int LENGTH_5K = 5000;
    private static final int LENGTH_100K = 100_000;

    public void addRead(final PrepRead read)
    {
        addToCounts(
                read.Chromosome, read.MateChromosome, read.start(), read.MatePosStart, read.orientation(), read.mateOrientation(),
                inferredInsertSizeAbs(read.record()));
    }

    public void processSupplementaryInfo(final ReadGroup readGroup)
    {
        // infer a short DEL or other category from the suppementary and its primary
        PrepRead read = null; // readGroup.reads().stream().filter(x -> x.hasSuppAlignment()).findFirst().orElse(null);

        for(PrepRead prepRead : readGroup.reads())
        {
            if(!prepRead.hasSuppAlignment())
                continue;

            if(read == null || read.isSupplementaryAlignment() && !prepRead.isSupplementaryAlignment())
                read = prepRead;
        }

        if(read == null)
            return;

        SupplementaryReadData suppData = read.supplementaryAlignment();

        int inferredDelLength = 0;
        Orientation suppOrientation = Orientation.fromByte(suppData.orientation());

        if(read.Chromosome.equals(suppData.Chromosome))
        {
            // only consider short DELs
            if(abs(suppData.Position - read.start()) > LENGTH_SHORT * 2)
                return;

            if(positionWithin(suppData.Position, read.start(), read.end()))
                return;

            boolean suppReadIsRightClipped;

            if(read.isSupplementaryAlignment())
                suppReadIsRightClipped = read.isRightClipped();
            else
                suppReadIsRightClipped = suppData.Cigar.endsWith("S");

            if(suppReadIsRightClipped)
            {
                // supp end position needs to be the inner part of a DEL
                if(!read.isLeftClipped() || read.leftClipLength() < read.rightClipLength())
                    return;

                int suppPosEnd = SamRecordUtils.getMateAlignmentEnd(suppData.Position, suppData.Cigar);

                if(suppData.Position >= read.start())
                    return;

                inferredDelLength = abs(read.start() - suppPosEnd);
            }
            else
            {
                if(!read.isRightClipped() || read.rightClipLength() < read.leftClipLength())
                    return;

                if(suppData.Position < read.end())
                    return;

                inferredDelLength = abs(suppData.Position - read.end());
            }

            if(inferredDelLength < 100 || inferredDelLength > LENGTH_SHORT)
                return;

            suppOrientation = read.orientation().opposite(); // to register as a DEL
        }

        addToCounts(
                read.Chromosome, suppData.Chromosome, read.start(), suppData.Position, read.orientation(),
                suppOrientation, inferredDelLength);
    }

    private void addToCounts(
            final String chromosome1, final String chromosome2, int posStart1, int posStart2,
            final Orientation orientation1, final Orientation orientation2, int insertSize)
    {
        if(!chromosome1.equals(chromosome2))
        {
            ++TypeCounts[DiscordantType.Translocation.ordinal()];
            return;
        }

        int distance = insertSize;

        if(orientation1 == orientation2)
        {
            if(distance < LENGTH_SHORT)
                TypeCounts[DiscordantType.InvLt1K.ordinal()] += 2; // assumes mate is in the same read group so not also registered
            else if(distance < LENGTH_5K)
                ++TypeCounts[DiscordantType.Inv1To5K.ordinal()];
            else if(distance < LENGTH_100K)
                ++TypeCounts[DiscordantType.Inv5To100K.ordinal()];
            else
                ++TypeCounts[DiscordantType.InvGt100K.ordinal()];
        }
        else
        {
            if((posStart1 < posStart2) == orientation1.isForward())
            {
                // DEL-type orientation
                if(distance < LENGTH_SHORT)
                    ++TypeCounts[DiscordantType.Del1Lt1K.ordinal()];
                else if(distance < LENGTH_5K)
                    ++TypeCounts[DiscordantType.Del1To5K.ordinal()];
                else if(distance < LENGTH_100K)
                    ++TypeCounts[DiscordantType.Del5To100K.ordinal()];
                else
                    ++TypeCounts[DiscordantType.DelGt100K.ordinal()];
            }
            else
            {
                if(distance < LENGTH_5K)
                    ++TypeCounts[DiscordantType.Dup1To5K.ordinal()];
                else if(distance < LENGTH_100K)
                    ++TypeCounts[DiscordantType.Dup5To100K.ordinal()];
                else
                    ++TypeCounts[DiscordantType.DupGt100K.ordinal()];
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
        long shortInvCount = TypeCounts[DiscordantType.InvLt1K.ordinal()] + TypeCounts[DiscordantType.Inv1To5K.ordinal()];
        return ProcessedReads > 0 ? shortInvCount / (double)ProcessedReads : 0;
    }

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

        return format("totalReads=%d writtenReads=%d %s", ProcessedReads, PrepReads, sj);
    }

    private static final String FLD_TOTAL_READS = "TotalReads";
    private static final String FLD_PREP_READS = "PrepReads";

    public static void writeDiscordantStats(
            final PrepConfig config, final long processedReads, final long writtenReads, final DiscordantStats discordantStats)
    {
        try
        {
            final String outputFileName = config.formFilename(DISCORDANT_STATS);

            BufferedWriter writer = createBufferedWriter(outputFileName, false);

            StringJoiner sj = new StringJoiner(TSV_DELIM);
            sj.add(FLD_TOTAL_READS);
            sj.add(FLD_PREP_READS);

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

    public static DiscordantStats loadDiscordantStats(final String filename)
    {
        try
        {
            List<String> lines = Files.readAllLines(Paths.get(filename));
            String[] values = lines.get(1).split(TSV_DELIM);

            long processedReads = Long.parseLong(values[0]);
            long prepReads = 0;

            long[] typeCounts = new long[DiscordantType.values().length];

            prepReads = Long.parseLong(values[1]);

            for(int i = 0; i < typeCounts.length; ++i)
            {
                typeCounts[i] = Long.parseLong(values[i + 2]);
            }

            return new DiscordantStats(processedReads, prepReads, typeCounts);
        }
        catch(Exception e)
        {
            SV_LOGGER.warn("failed to read discordant read statistics file({}): {}", filename, e.toString());
            return new DiscordantStats();
        }
    }
}
