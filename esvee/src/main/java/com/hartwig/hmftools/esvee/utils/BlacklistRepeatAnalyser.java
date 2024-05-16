package com.hartwig.hmftools.esvee.utils;

import static java.lang.Math.max;
import static java.lang.Math.min;
import static java.lang.Math.round;
import static java.lang.String.format;

import static com.hartwig.hmftools.common.genome.refgenome.RefGenomeSource.REF_GENOME;
import static com.hartwig.hmftools.common.genome.refgenome.RefGenomeSource.addRefGenomeConfig;
import static com.hartwig.hmftools.common.genome.refgenome.RefGenomeSource.loadRefGenome;
import static com.hartwig.hmftools.common.gripss.RepeatMaskAnnotations.REPEAT_MASK_FILE;
import static com.hartwig.hmftools.common.utils.config.ConfigUtils.addLoggingOptions;
import static com.hartwig.hmftools.common.utils.file.FileWriterUtils.addOutputOptions;
import static com.hartwig.hmftools.common.utils.file.FileWriterUtils.closeBufferedWriter;
import static com.hartwig.hmftools.common.utils.file.FileWriterUtils.createBufferedWriter;
import static com.hartwig.hmftools.esvee.AssemblyConfig.SV_LOGGER;
import static com.hartwig.hmftools.esvee.common.FileCommon.APP_NAME;
import static com.hartwig.hmftools.esvee.prep.PrepConfig.BLACKLIST_BED;
import static com.hartwig.hmftools.esvee.prep.PrepConfig.READ_LENGTH;
import static com.hartwig.hmftools.esvee.prep.PrepConstants.DEFAULT_READ_LENGTH;

import java.io.BufferedWriter;
import java.io.IOException;
import java.util.List;
import java.util.Set;

import com.google.common.collect.Lists;
import com.google.common.collect.Sets;
import com.hartwig.hmftools.common.codon.Nucleotides;
import com.hartwig.hmftools.common.genome.chromosome.HumanChromosome;
import com.hartwig.hmftools.common.genome.refgenome.RefGenomeInterface;
import com.hartwig.hmftools.common.genome.refgenome.RefGenomeVersion;
import com.hartwig.hmftools.common.gripss.RepeatMaskAnnotations;
import com.hartwig.hmftools.common.gripss.RepeatMaskData;
import com.hartwig.hmftools.common.region.BaseRegion;
import com.hartwig.hmftools.common.utils.config.ConfigBuilder;
import com.hartwig.hmftools.esvee.prep.BlacklistLocations;

import org.jetbrains.annotations.NotNull;

public class BlacklistRepeatAnalyser
{
    private final String mOutputFile;

    private final BlacklistLocations mBlacklistLocations;
    private final RepeatMaskAnnotations mRepeatMaskAnnotations;
    private final BufferedWriter mWriter;
    private final RefGenomeInterface mRefGenome;
    public final RefGenomeVersion mRefGenVersion;
    private final int mReadLength;
    private final double mMinRepeatPercent;

    private static final String MIN_REPEAT_PERC = "min_repeat_percent";
    private static final String OUTPUT_FILE = "output_file";

    private static final int MIN_REPEAT_BASES = 30;
    private static final double DEFAULT_MIN_REPEAT_PERC = 0.7;

    public BlacklistRepeatAnalyser(final ConfigBuilder configBuilder)
    {
        mBlacklistLocations = new BlacklistLocations(configBuilder.getValue(BLACKLIST_BED));
        mOutputFile = configBuilder.getValue(OUTPUT_FILE);
        mWriter = initialiseWriter();

        mRefGenome = loadRefGenome(configBuilder.getValue(REF_GENOME));
        mRefGenVersion = RefGenomeVersion.from(configBuilder);
        mReadLength = configBuilder.getInteger(READ_LENGTH);
        mMinRepeatPercent = configBuilder.getDecimal(MIN_REPEAT_PERC);

        mRepeatMaskAnnotations = new RepeatMaskAnnotations();
        if(configBuilder.hasValue(REPEAT_MASK_FILE))
        {
            if(!mRepeatMaskAnnotations.load(configBuilder.getValue(REPEAT_MASK_FILE), mRefGenVersion))
                System.exit(1);
        }

    }

    public void run()
    {
        SV_LOGGER.info("analysing repeats in {} blacklist locations", mBlacklistLocations.size());

        int regionCount = 0;
        for(HumanChromosome chromosome : HumanChromosome.values())
        {
            String chrStr = mRefGenVersion.versionedChromosome(chromosome.toString());

            List<BaseRegion> regions = mBlacklistLocations.getRegions(chrStr);

            if(regions == null)
                continue;

            for(BaseRegion region : regions)
            {
                processRegion(chrStr, region);
                ++regionCount;

                if((regionCount % 1000) == 0)
                {
                    SV_LOGGER.debug("processed {} regions", regionCount);
                }
            }
        }

        closeBufferedWriter(mWriter);

        SV_LOGGER.info("blacklist analysis complete");
    }

    private void processRegion(final String chromosome, final BaseRegion region)
    {
        List<RepeatData> repeats = Lists.newArrayList();

        String refBases = mRefGenome.getBaseString(chromosome, region.start(), region.end());
        int refBaseLength = refBases.length();
        int readShift = 1; // mReadLength / 10;

        RepeatData currentRepeat = null;
        int readIndex = 0;

        List<RepeatMaskData> rmMatches = mRepeatMaskAnnotations.findMatches(chromosome, region);

        while(readIndex < refBaseLength)
        {
            boolean atRegionEnd = readIndex + mReadLength >= refBaseLength;
            String readRefBases = refBases.substring(readIndex, min(readIndex + mReadLength, refBaseLength));

            if(readRefBases.length() < MIN_REPEAT_BASES)
                break;

            int[] baseCounts = calcNucleotideCounts(readRefBases);

            int minRepeatCount = max((int)round(mMinRepeatPercent * readRefBases.length()), MIN_REPEAT_BASES);

            char repeatedBase = 0;
            int repeatCount = 0;

            for(int b = 0; b < baseCounts.length; ++b)
            {
                if(baseCounts[b] >= minRepeatCount)
                {
                    repeatedBase = Nucleotides.DNA_BASES[b];
                    repeatCount = baseCounts[b];
                    break;
                }
            }

            if(repeatedBase > 0)
            {
                int repeatRegionStart = region.start() + readIndex;
                int repeatRegionEnd = repeatRegionStart + readRefBases.length() - 1;

                if(currentRepeat == null || repeatRegionStart > currentRepeat.Region.end())
                {
                    BaseRegion repeatRegion = new BaseRegion(repeatRegionStart, repeatRegionEnd);
                    currentRepeat = new RepeatData(repeatRegion, String.valueOf(repeatedBase));
                    currentRepeat.MinCount = currentRepeat.MaxCount = repeatCount;
                    repeats.add(currentRepeat);
                }
                else
                {
                    currentRepeat.Region.setEnd(repeatRegionEnd);
                    currentRepeat.MaxCount = max(currentRepeat.MaxCount, repeatCount);
                    currentRepeat.MinCount = min(currentRepeat.MinCount, repeatCount);
                }
            }

            if(atRegionEnd)
                break;

            readIndex += readShift;
        }

        if(!repeats.isEmpty() || !rmMatches.isEmpty())
            writeRegionData(chromosome, region, repeats, rmMatches);
    }

    private int[] calcNucleotideCounts(final String bases)
    {
        int[] baseCounts = new int[Nucleotides.DNA_BASES.length];

        for(int i = 0; i < bases.length(); ++i)
        {
            char base = bases.charAt(i);
            int baseIndex = Nucleotides.baseIndex(base);

            if(baseIndex >= 0)
                ++baseCounts[baseIndex];
        }

        return baseCounts;
    }

    private class RepeatData
    {
        public final BaseRegion Region;
        public String Repeat;
        public int MinCount;
        public int MaxCount;

        public RepeatData(final BaseRegion region, final String repeat)
        {
            Region = region;
            Repeat = repeat;
            MinCount = 0;
            MaxCount = 0;
        }
    }

    private void writeRegionData(
            final String chromosome, final BaseRegion region, final List<RepeatData> repeats, final List<RepeatMaskData> rmMatches)
    {
        try
        {
            String regionStr = format("%s\t%d\t%d", chromosome, region.start(), region.end());

            for(RepeatData repeatData : repeats)
            {
                mWriter.write(format("%s\t%s\t%s\t%d\t%d\t%d\t%s",
                        regionStr, "BASE_REPEAT", repeatData.Repeat, repeatData.Region.start(), repeatData.Region.end(),
                        repeatData.MaxCount, format("MinCount=%d",repeatData.MinCount)));
                mWriter.newLine();
            }

            // log once per class type
            Set<String> processedClasses = Sets.newHashSet();

            for(RepeatMaskData rmData : rmMatches)
            {
                if(processedClasses.contains(rmData.ClassType))
                    continue;

                processedClasses.add(rmData.ClassType);

                mWriter.write(format("%s\t%s\t%s\t%d\t%d\t%d\t%s",
                        regionStr, "REPEAT_MASK", rmData.ClassType,
                        rmData.Region.start(), rmData.Region.end(), 0, format("Id=%s", rmData.Id)));
                mWriter.newLine();
            }
        }
        catch(IOException e)
        {
            SV_LOGGER.error("failed to write output file: {}", e.toString());
        }
    }

    private BufferedWriter initialiseWriter()
    {
        try
        {
            BufferedWriter writer = createBufferedWriter(mOutputFile, false);

            writer.write("Chromosome\tPosStart\tPosEnd\tRepeatType\tRepeatInfo\tRepeatPosStart\tRepeatPosEnd\tCount\tOtherInfo");
            writer.newLine();

            return writer;
        }
        catch(IOException e)
        {
            SV_LOGGER.error("failed to initialise output file: {}", e.toString());
            return null;
        }
    }

    public static void main(@NotNull final String[] args)
    {
        ConfigBuilder configBuilder = new ConfigBuilder(APP_NAME);

        addRefGenomeConfig(configBuilder, true);
        configBuilder.addPath(BLACKLIST_BED, false, "Blacklist BED file");
        configBuilder.addConfigItem(OUTPUT_FILE, true, "Output filename");
        configBuilder.addInteger(READ_LENGTH, "Read length", DEFAULT_READ_LENGTH);
        configBuilder.addDecimal(MIN_REPEAT_PERC, "Min repeat perecent of read", DEFAULT_MIN_REPEAT_PERC);
        RepeatMaskAnnotations.addConfig(configBuilder);

        addOutputOptions(configBuilder);
        addLoggingOptions(configBuilder);

        configBuilder.checkAndParseCommandLine(args);

        BlacklistRepeatAnalyser blacklistRepeatAnalyser = new BlacklistRepeatAnalyser(configBuilder);
        blacklistRepeatAnalyser.run();
    }
}
