package com.hartwig.hmftools.isofox.fusion.cohort;

import static com.hartwig.hmftools.common.fusion.FusionCommon.FS_DOWN;
import static com.hartwig.hmftools.common.fusion.FusionCommon.FS_UP;
import static com.hartwig.hmftools.common.utils.ConfigUtils.addLoggingOptions;
import static com.hartwig.hmftools.common.utils.ConfigUtils.loadSampleIdsFile;
import static com.hartwig.hmftools.common.utils.ConfigUtils.setLogLevel;
import static com.hartwig.hmftools.common.utils.FileWriterUtils.OUTPUT_ID;
import static com.hartwig.hmftools.common.utils.FileWriterUtils.addOutputOptions;
import static com.hartwig.hmftools.common.utils.FileWriterUtils.closeBufferedWriter;
import static com.hartwig.hmftools.common.utils.FileWriterUtils.createBufferedWriter;
import static com.hartwig.hmftools.common.utils.FileWriterUtils.parseOutputDir;
import static com.hartwig.hmftools.common.utils.sv.StartEndIterator.SE_END;
import static com.hartwig.hmftools.common.utils.sv.StartEndIterator.SE_START;
import static com.hartwig.hmftools.isofox.IsofoxConfig.ISF_LOGGER;

import java.io.BufferedWriter;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.Arrays;
import java.util.List;
import java.util.Map;
import java.util.StringJoiner;

import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import com.hartwig.hmftools.isofox.fusion.FusionData;
import com.hartwig.hmftools.isofox.loader.DataLoaderConfig;

import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.CommandLineParser;
import org.apache.commons.cli.DefaultParser;
import org.apache.commons.cli.Options;
import org.apache.commons.cli.ParseException;
import org.jetbrains.annotations.NotNull;

public class FusionCompare
{
    private final String mFusionFileOrig;
    private final String mFusionFileNew;

    private final List<String> mSampleIds;
    private final BufferedWriter mWriter;

    private int mMatchedCount;
    private int mDiffCount;
    private int mUnmatchOrig;
    private int mUnmatchNew;

    private static final String FUSION_FILE_ORIG = "fusions_orig";
    private static final String FUSION_FILE_NEW = "fusions_new";
    private static final String SAMPLE_IDS_FILE = "sample_id_file";

    public FusionCompare(final CommandLine cmd)
    {
        mFusionFileOrig = cmd.getOptionValue(FUSION_FILE_ORIG);
        mFusionFileNew = cmd.getOptionValue(FUSION_FILE_NEW);
        mMatchedCount = 0;
        mDiffCount = 0;
        mUnmatchOrig = 0;
        mUnmatchNew = 0;

        mSampleIds = Lists.newArrayList();

        if(cmd.hasOption(SAMPLE_IDS_FILE))
        {
            mSampleIds.addAll(loadSampleIdsFile(cmd.getOptionValue(SAMPLE_IDS_FILE)));
        }

        String outputDir = parseOutputDir(cmd);
        mWriter = initialiseWriter(outputDir, cmd.getOptionValue(OUTPUT_ID));
    }

    public void run()
    {
        if(mWriter == null)
        {
            System.exit(1);
            return;
        }

        if(mSampleIds.isEmpty())
        {
            processSample(null, mFusionFileOrig, mFusionFileNew);
        }
        else
        {
            ISF_LOGGER.info("running comparisons for {} samples", mSampleIds.size());

            int processed = 0;
            for(String sampleId : mSampleIds)
            {
                String origFile = mFusionFileOrig.replaceAll("\\*", sampleId);
                String newFile = mFusionFileNew.replaceAll("\\*", sampleId);
                processSample(sampleId, origFile, newFile);

                ++processed;

                if(processed > 0 && (processed % 100) == 0)
                {
                    ISF_LOGGER.info("processed {} samples", processed);
                }
            }

            ISF_LOGGER.info("matched({}) diffs({}) unmatched(orig={} new={})",
                    mMatchedCount, mDiffCount, mUnmatchOrig, mUnmatchNew);
        }

        closeBufferedWriter(mWriter);
    }

    private void processSample(final String sampleId, final String origFile, final String newFile)
    {
        if(!Files.exists(Paths.get(origFile)))
        {
            ISF_LOGGER.error("invalid original fusions file({})", origFile);
            return;
        }

        if(!Files.exists(Paths.get(newFile)))
        {
            ISF_LOGGER.error("invalid new fusions file({})", newFile);
            return;
        }

        List<FusionData> origFusionsAll = FusionData.loadFromFile(Paths.get(origFile));
        List<FusionData> newFusionsAll = FusionData.loadFromFile(Paths.get(newFile));

        if(origFusionsAll.isEmpty() || newFusionsAll.isEmpty())
        {
            ISF_LOGGER.warn("empty fusion files: orig({}) new({})", origFusionsAll.size(), newFusionsAll.size());
            return;
        }

        if(sampleId != null)
        {
            ISF_LOGGER.debug("sample({}) loaded fusions orig({}) new({})", sampleId, origFusionsAll.size(), newFusionsAll.size());
        }
        else
        {
            ISF_LOGGER.info("loaded fusions orig({}) new({})", origFusionsAll.size(), newFusionsAll.size());
        }

        Map<String,List<FusionData>> origChrMap = buildChromosomeMap(origFusionsAll);
        Map<String,List<FusionData>> newChrMap = buildChromosomeMap(newFusionsAll);

        int matchedCount = 0;
        int diffCount = 0;
        int compared = 0;

        for(Map.Entry<String,List<FusionData>> entry : origChrMap.entrySet())
        {
            String chrPair = entry.getKey();
            List<FusionData> origFusions = entry.getValue();
            List<FusionData> newFusions = newChrMap.get(chrPair);

            if(newFusions == null)
                continue;

            int origIndex = 0;
            while(origIndex < origFusions.size())
            {
                FusionData origFusion = origFusions.get(origIndex);

                boolean matched = false;
                for(int newIndex = 0; newIndex < newFusions.size(); ++newIndex)
                {
                    FusionData newFusion = newFusions.get(newIndex);

                    if(matchOnCoords(origFusion, newFusion))
                    {
                        matched = true;
                        newFusions.remove(newIndex);

                        final List<String> diffs = Lists.newArrayList();
                        compareFusionDetails(origFusion, newFusion, diffs);

                        if(diffs.isEmpty())
                        {
                            ++matchedCount;
                        }
                        else
                        {
                            ++diffCount;
                            writeFusionDiffs(sampleId, "DIFFS", origFusion, newFusion, diffs);
                        }

                        break;
                    }
                }

                if(matched)
                    origFusions.remove(origIndex);
                else
                    ++origIndex;

                ++compared;

                if((compared % 1000) == 0)
                {
                    ISF_LOGGER.debug("compared {} fusions {}",
                            compared, sampleId != null ? String.format("for sample(%s)", sampleId) : "");
                }

                if(newFusions.isEmpty())
                    break;
            }
        }

        List<String> emptyDiffs = Lists.newArrayList();

        int unmatchOrig = 0;
        int unmatchNew = 0;

        for(List<FusionData> fusions : origChrMap.values())
        {
            for(FusionData fusion : fusions)
            {
                boolean isImmunePair = (fusion.GeneNames[FS_UP].startsWith("HLA") && fusion.GeneNames[FS_DOWN].startsWith("HLA"))
                        || (fusion.GeneNames[FS_UP].startsWith("IG") && fusion.GeneNames[FS_DOWN].startsWith("IG"));

                if(isImmunePair)
                {
                    writeFusionDiffs(sampleId, "IMMUNE_PAIR", fusion, null, emptyDiffs);
                }
                else
                {
                    writeFusionDiffs(sampleId, "NO_NEW", fusion, null, emptyDiffs);
                    ++unmatchOrig;
                }
            }
        }

        for(List<FusionData> fusions : newChrMap.values())
        {
            fusions.forEach(x -> writeFusionDiffs(sampleId, "NO_ORIG", null, x, emptyDiffs));
            unmatchNew += fusions.size();
        }

        if(sampleId != null)
        {
            ISF_LOGGER.debug("sample({}) matched({}) diffs({}) unmatched(orig={} new={})",
                    sampleId, matchedCount, diffCount, unmatchOrig, unmatchNew);
        }

        mMatchedCount += matchedCount;
        mUnmatchNew += unmatchNew;
        mUnmatchOrig += unmatchOrig;
        mDiffCount += diffCount;
    }

    private Map<String,List<FusionData>> buildChromosomeMap(final List<FusionData> fusions)
    {
        Map<String,List<FusionData>> chrMap = Maps.newHashMap();

        for(FusionData fusion : fusions)
        {
            String chrPair = fusion.Chromosomes[SE_START] + "_" + fusion.Chromosomes[SE_END];

            List<FusionData> chrFusions = chrMap.get(chrPair);
            if(chrFusions == null)
            {
                chrFusions = Lists.newArrayList();
                chrMap.put(chrPair, chrFusions);
            }

            chrFusions.add(fusion);
        }

        return chrMap;
    }

    private static boolean matchOnCoords(final FusionData fusion1, final FusionData fusion2)
    {
        if(!Arrays.equals(fusion1.Chromosomes, fusion2.Chromosomes))
            return false;

        if(!Arrays.equals(fusion1.JunctionPositions, fusion2.JunctionPositions))
            return false;

        if(!Arrays.equals(fusion1.JunctionOrientations, fusion2.JunctionOrientations))
            return false;

        return true;
    }

    private static void compareFusionDetails(final FusionData fusion1, final FusionData fusion2, final List<String> diffs)
    {
        if(fusion1.TotalFrags != fusion2.TotalFrags)
            diffs.add(String.format("totalFrags(%d/%d)", fusion1.TotalFrags, fusion2.TotalFrags));

        if(!Arrays.equals(fusion1.GeneNames, fusion2.GeneNames))
            diffs.add(String.format("genes(%s/%s)", fusion1.name(), fusion2.name()));
    }

    public BufferedWriter initialiseWriter(final String outputDir, final String outputId)
    {
        try
        {
            String outputFileName = outputDir + "isf.fusions_compare." + outputId + ".csv";

            final BufferedWriter writer = createBufferedWriter(outputFileName, false);

            if(!mSampleIds.isEmpty())
                writer.write("SampleId,");

            writer.write("MatchType,Name,IdOrig,IdNew,ChrUp,ChrDown,PosUp,PosDown,TotalFragsOrig,TotalFragsNew,SplitFragsOrig,SplitFragsNew,DiffInfo");
            writer.newLine();

            return writer;
        }
        catch(IOException e)
        {
            ISF_LOGGER.error("failed to initialise fusion comparison file: {}", e.toString());
            return null;
        }
    }

    private void writeFusionDiffs(
            final String sampleId, final String matchType, final FusionData origFusion, final FusionData newFusion, final List<String> diffs)
    {
        try
        {
            final FusionData refFusion = origFusion != null ? origFusion : newFusion;
            StringJoiner sj = new StringJoiner(";");
            diffs.forEach(x -> sj.add(x));

            if(sampleId != null)
                mWriter.write(String.format("%s,", sampleId));

            mWriter.write(String.format("%s,%s,%d,%d,%s,%s,%d,%d,%d,%d,%d,%d,%s",
                    matchType, origFusion != null ? origFusion.name() : newFusion.name(),
                    origFusion != null ? origFusion.Id : -1, newFusion != null ? newFusion.Id : -1,
                    refFusion.Chromosomes[SE_START], refFusion.Chromosomes[SE_END],
                    refFusion.JunctionPositions[SE_START], refFusion.JunctionPositions[SE_END],
                    origFusion != null ? origFusion.TotalFrags : -1, newFusion != null ? newFusion.TotalFrags : -1,
                    origFusion != null ? origFusion.SplitFrags : -1, newFusion != null ? newFusion.SplitFrags : -1,
                    sj.toString()));

            mWriter.newLine();
        }
        catch(IOException e)
        {
            ISF_LOGGER.error("failed to write fusion comparison file: {}", e.toString());
        }
    }

    public static void main(@NotNull final String[] args) throws ParseException
    {
        final Options options = DataLoaderConfig.createCmdLineOptions();

        options.addOption(FUSION_FILE_ORIG, true, "Original fusions file");
        options.addOption(FUSION_FILE_NEW, true, "New fusions file");
        options.addOption(SAMPLE_IDS_FILE, true, "Sample IDs file");
        addLoggingOptions(options);
        addOutputOptions(options);

        final CommandLineParser parser = new DefaultParser();
        final CommandLine cmd = parser.parse(options, args);

        setLogLevel(cmd);

        FusionCompare fusionCompare = new FusionCompare(cmd);
        fusionCompare.run();
    }

}
