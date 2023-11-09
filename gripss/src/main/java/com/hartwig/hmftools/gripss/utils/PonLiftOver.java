package com.hartwig.hmftools.gripss.utils;

import static java.lang.String.format;

import static com.hartwig.hmftools.common.genome.refgenome.GenomeLiftoverCache.UNMAPPED_POSITION;
import static com.hartwig.hmftools.common.genome.refgenome.RefGenomeVersion.V37;
import static com.hartwig.hmftools.common.genome.refgenome.RefGenomeVersion.V38;
import static com.hartwig.hmftools.common.utils.file.FileWriterUtils.createBufferedWriter;
import static com.hartwig.hmftools.common.utils.sv.SvCommonUtils.NEG_ORIENT_ID;
import static com.hartwig.hmftools.common.utils.sv.SvCommonUtils.POS_ORIENT;
import static com.hartwig.hmftools.common.utils.sv.SvCommonUtils.POS_ORIENT_ID;
import static com.hartwig.hmftools.gripss.GripssConfig.GR_LOGGER;

import java.io.BufferedWriter;
import java.io.IOException;
import java.util.List;
import java.util.Map;

import com.hartwig.hmftools.common.genome.chromosome.HumanChromosome;
import com.hartwig.hmftools.common.genome.refgenome.GenomeLiftoverCache;
import com.hartwig.hmftools.common.genome.refgenome.RefGenomeVersion;
import com.hartwig.hmftools.common.utils.config.ConfigBuilder;
import com.hartwig.hmftools.common.utils.config.ConfigUtils;
import com.hartwig.hmftools.gripss.pon.PonCache;
import com.hartwig.hmftools.gripss.pon.PonSglRegion;
import com.hartwig.hmftools.gripss.pon.PonSvRegion;

import org.jetbrains.annotations.NotNull;

public class PonLiftOver
{
    private final PonCache mPonCache;
    private final String mOutputSvFile;
    private final String mOutputSglFile;
    private final RefGenomeVersion mDestinationVersion;
    private final RefGenomeVersion mSourceVersion;
    private final GenomeLiftoverCache mLiftoverCache;

    private static final String OUTPUT_SV_PON_FILE = "output_pon_sv_file";
    private static final String OUTPUT_SGL_PON_FILE = "output_pon_sgl_file";
    private static final String DEST_REF_GENOME_VERSION = "dest_ref_genome_version";

    public PonLiftOver(final ConfigBuilder configBuilder)
    {
        mPonCache = new PonCache(configBuilder);
        mOutputSvFile = configBuilder.getValue(OUTPUT_SV_PON_FILE);
        mOutputSglFile = configBuilder.getValue(OUTPUT_SGL_PON_FILE);
        mDestinationVersion = RefGenomeVersion.from(configBuilder.getValue(DEST_REF_GENOME_VERSION));
        mSourceVersion = mDestinationVersion.is37() ? V38 : V37;

        mLiftoverCache = new GenomeLiftoverCache(true, mDestinationVersion == V38);
    }

    public void run()
    {
        if(!mPonCache.hasValidData() || mOutputSglFile == null || mOutputSvFile == null)
        {
            GR_LOGGER.error("invalid inputs");
            System.exit(1);
        }

        GR_LOGGER.info("Gripss PON lift-over");

        writeSvPonFile();
        writeSglPonFile();

        GR_LOGGER.info("Gripss PON lift-over complete");
    }

    private static final int LOG_COUNT = 1000000;

    private void writeSvPonFile()
    {
        GR_LOGGER.info("lifting over SV PON file to: {}", mOutputSvFile);

        int convertedRegions = 0;
        int droppedRegions = 0;

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
                    int posStartStart = mLiftoverCache.convertPosition(chrSourceStart, region.RegionStart.start());
                    int posStartEnd = mLiftoverCache.convertPosition(chrSourceStart, region.RegionStart.end());

                    String chrDestEnd = mDestinationVersion.versionedChromosome(region.RegionEnd.Chromosome);
                    int posEndStart =  mLiftoverCache.convertPosition(region.RegionEnd.Chromosome, region.RegionEnd.start());
                    int posEndEnd =  mLiftoverCache.convertPosition(region.RegionEnd.Chromosome, region.RegionEnd.end());

                    if(posStartStart == UNMAPPED_POSITION || posStartEnd == UNMAPPED_POSITION
                    || posEndStart == UNMAPPED_POSITION || posEndEnd == UNMAPPED_POSITION)
                    {
                        ++droppedRegions;
                        continue;
                    }

                    // fields: ChrStart,PosStartBegin,PosStartEnd,ChrEnd,PosEndBegin,PosEndEnd,Unknown,PonCount,OrientStart,OrientEnd
                    writer.write(format("%s\t%d\t%d\t%s\t%d\t%d\t%s\t%d\t%s\t%s",
                            chrDestStart, posStartStart, posStartEnd, chrDestEnd, posEndStart, posEndEnd, ".",
                            region.PonCount, orientToId(region.OrientStart), orientToId(region.OrientEnd)));

                    writer.newLine();

                    ++convertedRegions;

                    if((convertedRegions % LOG_COUNT) == 0)
                    {
                        GR_LOGGER.debug("converted {} positions", convertedRegions);
                    }
                }
            }

            writer.close();
        }
        catch(IOException e)
        {
            GR_LOGGER.error("failed to write SV PON output file: {}", e.toString());
            System.exit(1);
        }

        GR_LOGGER.info("coverted SV {} entries, dropped {}", convertedRegions, droppedRegions);
    }

    private void writeSglPonFile()
    {
        GR_LOGGER.info("lifting over SGL PON file to: {}", mOutputSvFile);

        int convertedRegions = 0;
        int droppedRegions = 0;

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
                    int posStartStart = mLiftoverCache.convertPosition(chrSourceStart, region.Region.start());
                    int posStartEnd = mLiftoverCache.convertPosition(chrSourceStart, region.Region.end());

                    if(posStartStart == UNMAPPED_POSITION || posStartEnd == UNMAPPED_POSITION)
                    {
                        ++droppedRegions;
                        continue;
                    }

                    // fields: Chr,PosBegin,PosEnd,Unknown,PonCount,Orientation
                    writer.write(format("%s\t%d\t%d\t%s\t%d\t%s",
                            chrDestStart, posStartStart, posStartEnd, ".", region.PonCount, orientToId(region.Orient)));

                    writer.newLine();

                    ++convertedRegions;

                    if((convertedRegions % LOG_COUNT) == 0)
                    {
                        GR_LOGGER.debug("converted {} positions", convertedRegions);
                    }
                }
            }

            writer.close();
        }
        catch(IOException e)
        {
            GR_LOGGER.error("failed to write SGL PON output file: {}", e.toString());
            System.exit(1);
        }

        GR_LOGGER.info("coverted SGL {} entries, dropped {}", convertedRegions, droppedRegions);
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
        PonCache.addConfig(configBuilder);

        ConfigUtils.addLoggingOptions(configBuilder);

        configBuilder.checkAndParseCommandLine(args);

        PonLiftOver ponLiftOver = new PonLiftOver(configBuilder);
        ponLiftOver.run();
    }
}
