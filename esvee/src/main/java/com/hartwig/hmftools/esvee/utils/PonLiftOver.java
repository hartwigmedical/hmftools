package com.hartwig.hmftools.esvee.utils;

import static com.hartwig.hmftools.common.genome.refgenome.GenomeLiftoverCache.UNMAPPED_POSITION;
import static com.hartwig.hmftools.common.genome.refgenome.RefGenomeVersion.V37;
import static com.hartwig.hmftools.common.genome.refgenome.RefGenomeVersion.V38;
import static com.hartwig.hmftools.common.utils.file.FileWriterUtils.createBufferedWriter;
import static com.hartwig.hmftools.common.utils.sv.SvCommonUtils.NEG_ORIENT_ID;
import static com.hartwig.hmftools.common.utils.sv.SvCommonUtils.POS_ORIENT;
import static com.hartwig.hmftools.common.utils.sv.SvCommonUtils.POS_ORIENT_ID;
import static com.hartwig.hmftools.common.utils.sv.SvCommonUtils.flipOrientation;
import static com.hartwig.hmftools.esvee.AssemblyConfig.SV_LOGGER;

import java.io.BufferedWriter;
import java.io.IOException;
import java.util.List;
import java.util.Map;

import com.hartwig.hmftools.common.genome.chromosome.HumanChromosome;
import com.hartwig.hmftools.common.genome.refgenome.GenomeLiftoverCache;
import com.hartwig.hmftools.common.genome.refgenome.RefGenomeVersion;
import com.hartwig.hmftools.common.utils.config.ConfigBuilder;
import com.hartwig.hmftools.common.utils.config.ConfigUtils;
import com.hartwig.hmftools.esvee.caller.annotation.PonCache;
import com.hartwig.hmftools.esvee.caller.annotation.PonSvRegion;
import com.hartwig.hmftools.esvee.caller.annotation.PonSglRegion;

import org.jetbrains.annotations.NotNull;

public class PonLiftOver
{
    private final PonCache mPonCache;
    private final String mOutputSvFile;
    private final String mOutputSglFile;
    private final RefGenomeVersion mDestinationVersion;
    private final RefGenomeVersion mSourceVersion;
    private final GenomeLiftoverCache mLiftoverCache;
    private final int mMaxLengthDiff;

    private int mConverted;
    private int mFailedLiftover;
    private int mFailedMaxLengthDiff;

    private static final String OUTPUT_SV_PON_FILE = "output_pon_sv_file";
    private static final String OUTPUT_SGL_PON_FILE = "output_pon_sgl_file";
    private static final String DEST_REF_GENOME_VERSION = "dest_ref_genome_version";
    private static final String MAX_LENGTH_DIFF = "max_length_diff";

    public PonLiftOver(final ConfigBuilder configBuilder)
    {
        mPonCache = new PonCache(configBuilder);
        mOutputSvFile = configBuilder.getValue(OUTPUT_SV_PON_FILE);
        mOutputSglFile = configBuilder.getValue(OUTPUT_SGL_PON_FILE);
        mDestinationVersion = RefGenomeVersion.from(configBuilder.getValue(DEST_REF_GENOME_VERSION));
        mSourceVersion = mDestinationVersion.is37() ? V38 : V37;
        mMaxLengthDiff = configBuilder.getInteger(MAX_LENGTH_DIFF);

        mLiftoverCache = new GenomeLiftoverCache(true);

        mConverted = 0;
        mFailedLiftover = 0;
        mFailedMaxLengthDiff = 0;
    }

    public void run()
    {
        if(!mPonCache.hasValidData() || mOutputSglFile == null || mOutputSvFile == null)
        {
            SV_LOGGER.error("invalid inputs");
            System.exit(1);
        }

        SV_LOGGER.info("Gripss PON lift-over");

        writeSvPonFile();
        writeSglPonFile();

        SV_LOGGER.info("converted {} entries, failed liftover({}) maxLengthDiff({})",
                mConverted, mFailedLiftover, mFailedMaxLengthDiff);

        SV_LOGGER.info("Gripss PON lift-over complete");
    }

    private static final int LOG_COUNT = 1000000;

    private void writeSvPonFile()
    {
        SV_LOGGER.info("lifting over SV PON file to: {}", mOutputSvFile);

        try
        {
            BufferedWriter writer = createBufferedWriter(mOutputSvFile, false);

            Map<String,List<PonSvRegion>> chrSvRegionsMap = mPonCache.svRegions();

            for(HumanChromosome chromosome : HumanChromosome.values())
            {
                String chrSourceStart = mSourceVersion.versionedChromosome(chromosome.toString());
                String chrDestStart = mDestinationVersion.versionedChromosome(chrSourceStart);

                List<PonSvRegion> regions = chrSvRegionsMap.get(chrSourceStart);

                if(regions == null)
                    continue;

                for(PonSvRegion region : regions)
                {
                    int posStartStart = mLiftoverCache.convertPosition(chrSourceStart, region.RegionStart.start(), mDestinationVersion);
                    int posStartEnd = mLiftoverCache.convertPosition(chrSourceStart, region.RegionStart.end(), mDestinationVersion);

                    String chrDestEnd = mDestinationVersion.versionedChromosome(region.RegionEnd.Chromosome);
                    int posEndStart =  mLiftoverCache.convertPosition(region.RegionEnd.Chromosome, region.RegionEnd.start(), mDestinationVersion);
                    int posEndEnd =  mLiftoverCache.convertPosition(region.RegionEnd.Chromosome, region.RegionEnd.end(), mDestinationVersion);

                    if(posStartStart == UNMAPPED_POSITION || posStartEnd == UNMAPPED_POSITION
                    || posEndStart == UNMAPPED_POSITION || posEndEnd == UNMAPPED_POSITION)
                    {
                        ++mFailedLiftover;
                        continue;
                    }

                    byte orientStart = region.OrientStart.asByte();
                    byte orientEnd = region.OrientEnd.asByte();

                    // check for changes to direction for SVs
                    if(posStartStart > posStartEnd)
                    {
                        int tmp = posStartEnd;
                        posStartEnd = posStartStart;
                        posStartStart = tmp;
                        orientStart = flipOrientation(orientStart);
                    }

                    if(posEndStart > posEndEnd)
                    {
                        int tmp = posEndEnd;
                        posEndEnd = posEndStart;
                        posEndStart = tmp;
                        orientEnd = flipOrientation(orientEnd);
                    }

                    // check for substantial changes in region length, indicating a possible error in lift-over
                    int lengthStartNew = posStartEnd - posStartStart;
                    int lengthEndNew = posEndEnd - posEndStart;

                    if(lengthStartNew - region.RegionStart.length() > mMaxLengthDiff
                    || lengthEndNew - region.RegionEnd.length() > mMaxLengthDiff)
                    {
                        /*
                        GR_LOGGER.debug("LENGTH: SV_START,{},{},{},{},{},{},{}",
                                chrDestStart, region.RegionStart.start(), region.RegionStart.end(), region.RegionStart.length(),
                                posStartStart, posStartEnd, posStartEnd - posStartStart);
                        */
                        ++mFailedMaxLengthDiff;
                        continue;

                    }

                    if(chrDestStart.equals(chrDestEnd) && posStartStart > posEndStart)
                    {
                        SV_LOGGER.trace("swapping start region({}:{}-{}) with end region({}:{}-{})",
                                chrDestStart, posStartStart, posStartEnd, chrDestEnd, posEndStart, posEndEnd);

                        // swap coords if the start region is now after the end region
                        writer.write(String.format("%s\t%d\t%d\t%s\t%d\t%d\t%s\t%d\t%s\t%s",
                                chrDestStart, posEndStart, posEndEnd, chrDestEnd, posStartStart, posStartEnd, ".",
                                region.PonCount, orientToId(orientEnd), orientToId(orientStart)));
                    }
                    else
                    {
                        // fields: ChrStart,PosStartBegin,PosStartEnd,ChrEnd,PosEndBegin,PosEndEnd,Unknown,PonCount,OrientStart,OrientEnd
                        writer.write(String.format("%s\t%d\t%d\t%s\t%d\t%d\t%s\t%d\t%s\t%s",
                                chrDestStart, posStartStart, posStartEnd, chrDestEnd, posEndStart, posEndEnd, ".",
                                region.PonCount, orientToId(orientStart), orientToId(orientEnd)));
                    }

                    writer.newLine();

                    ++mConverted;

                    if((mConverted % LOG_COUNT) == 0)
                    {
                        SV_LOGGER.debug("converted {} entries", mConverted);
                    }
                }
            }

            writer.close();
        }
        catch(IOException e)
        {
            SV_LOGGER.error("failed to write SV PON output file: {}", e.toString());
            System.exit(1);
        }
    }

    private void writeSglPonFile()
    {
        SV_LOGGER.info("lifting over SGL PON file to: {}", mOutputSvFile);

        try
        {
            BufferedWriter writer = createBufferedWriter(mOutputSglFile, false);

            Map<String,List<PonSglRegion>> chrSvRegionsMap = mPonCache.sglRegions();

            for(HumanChromosome chromosome : HumanChromosome.values())
            {
                String chrSourceStart = mSourceVersion.versionedChromosome(chromosome.toString());
                String chrDestStart = mDestinationVersion.versionedChromosome(chrSourceStart);

                List<PonSglRegion> regions = chrSvRegionsMap.get(chrSourceStart);

                if(regions == null)
                    continue;

                for(PonSglRegion region : regions)
                {
                    int posStartStart = mLiftoverCache.convertPosition(chrSourceStart, region.Region.start(), mDestinationVersion);
                    int posStartEnd = mLiftoverCache.convertPosition(chrSourceStart, region.Region.end(), mDestinationVersion);
                    byte orientStart = region.Orient.asByte();

                    if(posStartStart == UNMAPPED_POSITION || posStartEnd == UNMAPPED_POSITION)
                    {
                        ++mFailedLiftover;
                        continue;
                    }

                    if(posStartStart > posStartEnd)
                    {
                        int tmp = posStartEnd;
                        posStartEnd = posStartStart;
                        posStartStart = tmp;
                        orientStart = flipOrientation(orientStart);
                    }

                    int lengthNew = posStartEnd - posStartStart;

                    if(lengthNew - region.Region.length() > mMaxLengthDiff)
                    {
                        /*
                        GR_LOGGER.debug("LENGTH: SGL,{},{},{},{},{},{},{}",
                                chrDestStart, region.Region.start(), region.Region.end(), region.Region.length(),
                                posStartStart, posStartEnd, posStartEnd - posStartStart);
                        */

                        ++mFailedMaxLengthDiff;
                        continue;
                    }

                    // fields: Chr,PosBegin,PosEnd,Unknown,PonCount,Orientation
                    writer.write(String.format("%s\t%d\t%d\t%s\t%d\t%s",
                            chrDestStart, posStartStart, posStartEnd, ".", region.PonCount, orientToId(orientStart)));

                    writer.newLine();

                    ++mConverted;

                    if((mConverted % LOG_COUNT) == 0)
                    {
                        SV_LOGGER.debug("converted {} entries", mConverted);
                    }
                }
            }

            writer.close();
        }
        catch(IOException e)
        {
            SV_LOGGER.error("failed to write SGL PON output file: {}", e.toString());
            System.exit(1);
        }
    }

    private static String orientToId(final byte orientation)
    {
        return orientation == POS_ORIENT ? POS_ORIENT_ID : NEG_ORIENT_ID;
    }

    public static void main(@NotNull final String[] args)
    {
        ConfigBuilder configBuilder = new ConfigBuilder();
        configBuilder.addConfigItem(OUTPUT_SGL_PON_FILE, true, "Output SGL PON file");
        configBuilder.addConfigItem(OUTPUT_SV_PON_FILE, true, "Output SV PON file");
        configBuilder.addConfigItem(DEST_REF_GENOME_VERSION, true, "Destination ref genome version");
        configBuilder.addInteger(MAX_LENGTH_DIFF, "Permitted lifted region difference in length vs original", 5);
        PonCache.addConfig(configBuilder);

        ConfigUtils.addLoggingOptions(configBuilder);

        configBuilder.checkAndParseCommandLine(args);

        PonLiftOver ponLiftOver = new PonLiftOver(configBuilder);
        ponLiftOver.run();
    }
}
