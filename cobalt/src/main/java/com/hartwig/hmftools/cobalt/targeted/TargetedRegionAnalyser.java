package com.hartwig.hmftools.cobalt.targeted;

import static java.lang.Math.round;
import static java.lang.String.format;

import static com.hartwig.hmftools.cobalt.CobaltConfig.CB_LOGGER;
import static com.hartwig.hmftools.common.genome.bed.BedFileReader.loadBedFile;
import static com.hartwig.hmftools.common.genome.gc.GCProfileFactory.GC_PROFILE;
import static com.hartwig.hmftools.common.genome.gc.GCProfileFactory.addGcProfilePath;
import static com.hartwig.hmftools.common.genome.refgenome.RefGenomeSource.REF_GENOME;
import static com.hartwig.hmftools.common.genome.refgenome.RefGenomeSource.REF_GENOME_CFG_DESC;
import static com.hartwig.hmftools.common.genome.refgenome.RefGenomeSource.loadRefGenome;
import static com.hartwig.hmftools.common.utils.config.CommonConfig.TARGET_REGIONS_BED;
import static com.hartwig.hmftools.common.utils.config.CommonConfig.TARGET_REGIONS_BED_DESC;
import static com.hartwig.hmftools.common.utils.config.ConfigUtils.addLoggingOptions;
import static com.hartwig.hmftools.common.utils.config.ConfigUtils.setLogLevel;
import static com.hartwig.hmftools.common.utils.file.FileWriterUtils.addOutputDir;
import static com.hartwig.hmftools.common.utils.file.FileWriterUtils.createBufferedWriter;
import static com.hartwig.hmftools.common.utils.file.FileWriterUtils.parseOutputDir;

import java.io.BufferedWriter;
import java.io.IOException;
import java.util.List;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.cobalt.norm.GcProfileCache;
import com.hartwig.hmftools.common.genome.gc.GCProfile;
import com.hartwig.hmftools.common.genome.gc.GcCalcs;
import com.hartwig.hmftools.common.genome.refgenome.RefGenomeSource;
import com.hartwig.hmftools.common.utils.config.ConfigBuilder;
import com.hartwig.hmftools.common.utils.sv.ChrBaseRegion;

public class TargetedRegionAnalyser
{
    private final List<ChrBaseRegion> mTargetRegions;
    private final GcProfileCache mGcProfileCache;
    private final RefGenomeSource mRefGenome;
    private final String mOutputDir;

    public TargetedRegionAnalyser(final ConfigBuilder configBuilder)
    {
        mTargetRegions = Lists.newArrayList();

        try
        {
            mTargetRegions.addAll(loadBedFile(configBuilder.getValue(TARGET_REGIONS_BED)));
        }
        catch(Exception e)
        {
            System.exit(1);
        }

        mGcProfileCache = new GcProfileCache(configBuilder.getValue(GC_PROFILE));

        String refGenomeFile = configBuilder.getValue(REF_GENOME);
        mRefGenome = loadRefGenome(refGenomeFile);

        mOutputDir = parseOutputDir(configBuilder);
    }

    public void run()
    {
        CB_LOGGER.info("running Cobalt targeted region GC analyser");

        String outputFile = mOutputDir + "target_regions_gc_data.tsv";

        try
        {
            BufferedWriter writer = createBufferedWriter(outputFile, false);

            writer.write(format("Chromosome\tPosStart\tPosEnd\tBucketGc\tRegionGc\tBucketMappability"));
            writer.newLine();

            for(ChrBaseRegion region : mTargetRegions)
            {
                int bucketStart = round(region.start() / 1000) * 1000 + 1;
                GCProfile gcProfile = mGcProfileCache.findGcProfile(region.Chromosome, bucketStart);

                if(gcProfile == null)
                {
                    CB_LOGGER.error("region({}) bucket({}) GC profile entry not found", region, bucketStart);
                    continue;
                }

                String refSequence = mRefGenome.getBaseString(region.Chromosome, region.start(), region.end());

                double seqGcRatio = GcCalcs.calcGcPercent(refSequence);

                writer.write(format("%s\t%d\t%d\t%.3f\t%.3f\t%.2f",
                        region.Chromosome, region.start(), region.end(), gcProfile.gcContent(),
                        seqGcRatio, gcProfile.mappablePercentage()));

                writer.newLine();
            }

            writer.close();
        }
        catch(IOException e)
        {
            e.printStackTrace();
            CB_LOGGER.error("failed to write target regions GC data", e.toString());
            System.exit(1);
        }

        CB_LOGGER.info("Cobalt normalisation file generation complete");
    }

    public static void main(final String[] args)
    {
        ConfigBuilder configBuilder = new ConfigBuilder("Cobalt");

        configBuilder.addPath(TARGET_REGIONS_BED, true, TARGET_REGIONS_BED_DESC);
        configBuilder.addPath(REF_GENOME, false, REF_GENOME_CFG_DESC);
        addGcProfilePath(configBuilder, true);
        addOutputDir(configBuilder);
        addLoggingOptions(configBuilder);

        configBuilder.checkAndParseCommandLine(args);

        setLogLevel(configBuilder);

        TargetedRegionAnalyser targetedRegionAnalyser = new TargetedRegionAnalyser(configBuilder);
        targetedRegionAnalyser.run();
    }
}
