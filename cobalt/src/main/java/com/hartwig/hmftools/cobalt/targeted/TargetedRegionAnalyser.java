package com.hartwig.hmftools.cobalt.targeted;

import static java.lang.Math.max;
import static java.lang.Math.min;

import static com.hartwig.hmftools.cobalt.CobaltConfig.CB_LOGGER;
import static com.hartwig.hmftools.common.genome.gc.GCProfileFactory.GC_PROFILE;
import static com.hartwig.hmftools.common.genome.gc.GCProfileFactory.addGcProfilePath;
import static com.hartwig.hmftools.common.genome.refgenome.RefGenomeSource.REF_GENOME;
import static com.hartwig.hmftools.common.genome.refgenome.RefGenomeSource.REF_GENOME_CFG_DESC;
import static com.hartwig.hmftools.common.genome.refgenome.RefGenomeSource.loadRefGenome;
import static com.hartwig.hmftools.common.region.BedFileReader.loadBedFile;
import static com.hartwig.hmftools.common.utils.config.CommonConfig.TARGET_REGIONS_BED;
import static com.hartwig.hmftools.common.utils.config.CommonConfig.TARGET_REGIONS_BED_DESC;
import static com.hartwig.hmftools.common.utils.config.ConfigUtils.addLoggingOptions;
import static com.hartwig.hmftools.common.utils.config.ConfigUtils.setLogLevel;
import static com.hartwig.hmftools.common.utils.file.FileDelimiters.ITEM_DELIM;
import static com.hartwig.hmftools.common.utils.file.FileDelimiters.TSV_DELIM;
import static com.hartwig.hmftools.common.utils.file.FileWriterUtils.addOutputDir;
import static com.hartwig.hmftools.common.utils.file.FileWriterUtils.closeBufferedWriter;
import static com.hartwig.hmftools.common.utils.file.FileWriterUtils.createBufferedWriter;
import static com.hartwig.hmftools.common.utils.file.FileWriterUtils.parseOutputDir;

import java.io.BufferedWriter;
import java.io.IOException;
import java.util.List;
import java.util.Map;
import java.util.StringJoiner;
import java.util.stream.Collectors;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.cobalt.norm.GcProfileCache;
import com.hartwig.hmftools.common.genome.gc.GCProfile;
import com.hartwig.hmftools.common.genome.gc.GcCalcs;
import com.hartwig.hmftools.common.genome.refgenome.RefGenomeSource;
import com.hartwig.hmftools.common.region.BaseRegion;
import com.hartwig.hmftools.common.region.ChrBaseRegion;
import com.hartwig.hmftools.common.utils.config.ConfigBuilder;

public class TargetedRegionAnalyser
{
    private final List<ChrBaseRegion> mTargetRegions;
    private final GcProfileCache mGcProfileCache;
    private final RefGenomeSource mRefGenome;
    private final String mOutputDir;
    private final BufferedWriter mWriter;

    private static final int PROBE_LENGTH = 120;

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

        mWriter = initialiseWriter();
    }

    public void run()
    {
        CB_LOGGER.info("running Cobalt targeted region GC analyser");

        Map<String, List<GCProfile>> chrGcProfileMap = mGcProfileCache.chrGcProfiles();

        for(Map.Entry<String, List<GCProfile>> entry : chrGcProfileMap.entrySet())
        {
            String chromosome = entry.getKey();
            List<GCProfile> gcProfiles = entry.getValue();

            List<ChrBaseRegion> regions = mTargetRegions.stream().filter(x -> x.Chromosome.equals(chromosome)).toList();

            if(regions.isEmpty())
            {
                continue;
            }

            int regionIndex = 0;
            ChrBaseRegion currentRegion = regions.get(regionIndex);

            for(GCProfile gcProfile : gcProfiles)
            {
                if(gcProfile.end() < currentRegion.start())
                {
                    continue;
                }

                int regionBaseCount = 0;
                StringBuilder regionBaseSequence = new StringBuilder();
                List<BaseRegion> overlappingRegions = Lists.newArrayList();
                boolean regionOverlapsWindows = false;

                while(currentRegion.start() <= gcProfile.end())
                {
                    int regionOverlapStart = max(currentRegion.start(), gcProfile.start());
                    int regionOverlapEnd = min(currentRegion.end(), gcProfile.end());

                    int overlapCount = regionOverlapEnd - regionOverlapStart + 1;

                    if(overlapCount < PROBE_LENGTH)
                    {
                        if(currentRegion.baseLength() >= PROBE_LENGTH)
                        {
                            if(currentRegion.start() < gcProfile.start())
                            {
                                regionOverlapStart = regionOverlapEnd - PROBE_LENGTH + 1;
                            }
                            else
                            {
                                // region overlaps with next window
                                regionOverlapEnd = regionOverlapStart + PROBE_LENGTH - 1;
                            }
                        }
                        else
                        {
                            // expand from centre of region to the required probe length
                            int requiredExtraFlankingBases = (PROBE_LENGTH - currentRegion.baseLength()) / 2;
                            regionOverlapStart = currentRegion.start() - requiredExtraFlankingBases;
                            regionOverlapEnd = currentRegion.end() + requiredExtraFlankingBases;
                        }

                        overlapCount = regionOverlapEnd - regionOverlapStart + 1;
                    }

                    if(regionOverlapStart > regionOverlapEnd)
                    {
                        CB_LOGGER.error("invalid overlap");
                    }

                    try
                    {
                        String overlapSequence = mRefGenome.getBaseString(chromosome, regionOverlapStart, regionOverlapEnd);
                        regionBaseSequence.append(overlapSequence);
                    }
                    catch(Exception e)
                    {
                        CB_LOGGER.error("invalid ref sequence");
                    }

                    regionBaseCount += overlapCount;
                    overlappingRegions.add(new BaseRegion(regionOverlapStart, regionOverlapEnd));

                    regionOverlapsWindows = currentRegion.end() > gcProfile.end();

                    ++regionIndex;

                    if(regionIndex >= regions.size())
                    {
                        break;
                    }

                    currentRegion = regions.get(regionIndex);
                }

                double regionGc = GcCalcs.calcGcPercent(regionBaseSequence.toString());

                writeGcWindowData(gcProfile, regionBaseCount, regionGc, overlappingRegions);

                if(regionOverlapsWindows)
                {
                    --regionIndex;
                    currentRegion = regions.get(regionIndex);
                }

                if(regionIndex >= regions.size())
                {
                    break;
                }
            }
        }

        closeBufferedWriter(mWriter);

        CB_LOGGER.info("Cobalt normalisation file generation complete");
    }

    private BufferedWriter initialiseWriter()
    {
        try
        {
            String outputFile = mOutputDir + "target_regions_gc_data.tsv";
            BufferedWriter writer = createBufferedWriter(outputFile, false);

            writer.write("Chromosome\tGcWindowStart\tGcWindowEnd\tWindowGc\tRegionBaseCount\tRegionGc\tRegionIntervals\tWindowMappability");
            writer.newLine();
            return writer;
        }
        catch(IOException e)
        {
            CB_LOGGER.error("failed to write target regions GC data", e);
            System.exit(1);
        }

        return null;
    }

    private void writeGcWindowData(final GCProfile gcProfile, int regionBaseCount, double regionGc,
            final List<BaseRegion> overlappingRegions)
    {
        if(mWriter == null)
        {
            return;
        }

        try
        {
            StringJoiner sj = new StringJoiner(TSV_DELIM);
            sj.add(gcProfile.chromosome());
            sj.add(String.valueOf(gcProfile.start()));
            sj.add(String.valueOf(gcProfile.end()));
            sj.add(String.valueOf(gcProfile.gcContent()));
            sj.add(String.valueOf(regionBaseCount));
            sj.add(String.valueOf(regionGc));

            StringJoiner intervalsSj = new StringJoiner(ITEM_DELIM);
            overlappingRegions.forEach(x -> intervalsSj.add(x.toString()));
            sj.add(intervalsSj.toString());

            sj.add(String.valueOf(gcProfile.mappablePercentage()));
            mWriter.write(sj.toString());
            mWriter.newLine();
        }
        catch(IOException e)
        {
            CB_LOGGER.error("failed to write target regions GC data", e);
            System.exit(1);
        }
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
