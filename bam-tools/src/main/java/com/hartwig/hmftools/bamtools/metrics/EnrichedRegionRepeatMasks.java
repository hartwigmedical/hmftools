package com.hartwig.hmftools.bamtools.metrics;

import static java.lang.Math.max;
import static java.lang.Math.min;
import static java.lang.String.format;

import static com.hartwig.hmftools.bamtools.common.CommonUtils.APP_NAME;
import static com.hartwig.hmftools.bamtools.common.CommonUtils.BT_LOGGER;
import static com.hartwig.hmftools.bamtools.metrics.OffTargetFragments.FLD_PEAK_POS_END;
import static com.hartwig.hmftools.bamtools.metrics.OffTargetFragments.FLD_PEAK_POS_START;
import static com.hartwig.hmftools.common.genome.refgenome.RefGenomeSource.addRefGenomeVersion;
import static com.hartwig.hmftools.common.gripss.RepeatMaskAnnotations.REPEAT_MASK_FILE;
import static com.hartwig.hmftools.common.utils.TaskExecutor.addThreadOptions;
import static com.hartwig.hmftools.common.utils.TaskExecutor.parseThreads;
import static com.hartwig.hmftools.common.utils.config.ConfigUtils.addLoggingOptions;
import static com.hartwig.hmftools.common.utils.file.CommonFields.FLD_CHROMOSOME;
import static com.hartwig.hmftools.common.utils.file.FileDelimiters.TSV_DELIM;
import static com.hartwig.hmftools.common.utils.file.FileDelimiters.TSV_EXTENSION;
import static com.hartwig.hmftools.common.utils.file.FileReaderUtils.createFieldsIndexMap;
import static com.hartwig.hmftools.common.utils.file.FileWriterUtils.closeBufferedWriter;
import static com.hartwig.hmftools.common.utils.file.FileWriterUtils.createBufferedReader;
import static com.hartwig.hmftools.common.utils.file.FileWriterUtils.createBufferedWriter;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.List;
import java.util.Map;
import java.util.StringJoiner;
import java.util.concurrent.Callable;
import java.util.stream.Collectors;

import com.hartwig.hmftools.common.genome.refgenome.RefGenomeVersion;
import com.hartwig.hmftools.common.gripss.RepeatMaskAnnotations;
import com.hartwig.hmftools.common.gripss.RepeatMaskData;
import com.hartwig.hmftools.common.region.ChrBaseRegion;
import com.hartwig.hmftools.common.utils.TaskExecutor;
import com.hartwig.hmftools.common.utils.config.ConfigBuilder;

import org.jetbrains.annotations.NotNull;

public class EnrichedRegionRepeatMasks
{
    private final String mOutputFile;
    private final String mInputRegionsFile;
    private final RepeatMaskAnnotations mRmAnnotation;
    private BufferedWriter mWriter;
    private final int mThreads;

    private static final String INPUT_REGION_FILE = "input_region_file";
    private static final String OUTPUT_FILE = "output_file";

    public EnrichedRegionRepeatMasks(final ConfigBuilder configBuilder)
    {
        mRmAnnotation = new RepeatMaskAnnotations();

        RefGenomeVersion refGenomeVersion = RefGenomeVersion.from(configBuilder);

        if(!mRmAnnotation.load(configBuilder.getValue(REPEAT_MASK_FILE), refGenomeVersion))
            System.exit(1);

        mInputRegionsFile = configBuilder.getValue(INPUT_REGION_FILE);
        mOutputFile = configBuilder.getValue(OUTPUT_FILE);

        mWriter = null;
        mThreads = parseThreads(configBuilder);
    }

    public void run()
    {
        BT_LOGGER.info("running enriched region repeat-mask annotations");

        List<RegionData> regions = loadRegions();

        List<AnnotationTask> annotationTasks = new ArrayList<>();

        for(int i = 0; i < mThreads; ++i)
        {
            annotationTasks.add(new AnnotationTask(i));
        }

        int taskIndex = 0;
        for(RegionData region : regions)
        {
            if(taskIndex >= annotationTasks.size())
                taskIndex = 0;

            annotationTasks.get(taskIndex).Regions.add(region);

            ++taskIndex;
        }

        final List<Callable> callableList = annotationTasks.stream().collect(Collectors.toList());
        if(!TaskExecutor.executeTasks(callableList, mThreads))
            System.exit(1);

        closeBufferedWriter(mWriter);

        BT_LOGGER.info("Enriched regions repeat-masker annotation complete");
    }

    private class RegionData extends ChrBaseRegion
    {
        public final String Data;

        public RegionData(final String chromosome, final int posStart, final int posEnd, final String data
        )
        {
            super(chromosome, posStart, posEnd);
            Data = data;
        }
    }

    private List<RegionData> loadRegions()
    {
        List<RegionData> regions = new ArrayList<>();

        try
        {
            BufferedReader fileReader = createBufferedReader(mInputRegionsFile);

            String header = fileReader.readLine();

            String outputFile = mOutputFile != null ?
                    mOutputFile : mInputRegionsFile.replaceAll(TSV_EXTENSION, ".repeat_mask" + TSV_EXTENSION);

            mWriter = createBufferedWriter(outputFile, false);

            StringJoiner sj = new StringJoiner(TSV_DELIM);
            sj.add(header);
            sj.add("RmRegions").add("RmOverlapBases").add("RmClassTypes").add("RmRepeat");

            mWriter.write(sj.toString());
            mWriter.newLine();

            final Map<String,Integer> fieldsIndexMap = createFieldsIndexMap(header, TSV_DELIM);

            int chrIndex = fieldsIndexMap.get(FLD_CHROMOSOME);
            int posStartIndex = fieldsIndexMap.get(FLD_PEAK_POS_START);
            int posEndIndex = fieldsIndexMap.get(FLD_PEAK_POS_END);

            String line = null;

            while((line = fileReader.readLine()) != null)
            {
                final String[] values = line.split(TSV_DELIM, -1);

                RegionData regionData = new RegionData(
                        values[chrIndex], Integer.parseInt(values[posStartIndex]), Integer.parseInt(values[posEndIndex]), line);

                regions.add(regionData);
            }

            BT_LOGGER.info("loaded {} region entries from file({})", regions.size(), mInputRegionsFile);
        }
        catch(IOException e)
        {
            BT_LOGGER.error("failed to load regions from file({}): {}", mInputRegionsFile, e.toString());
            System.exit(1);
        }

        return regions;
    }

    private class AnnotationTask implements Callable
    {
        private final int mTaskId;
        public final List<RegionData> Regions;

        public AnnotationTask(int taskId)
        {
            mTaskId = taskId;
            Regions = new ArrayList<>();
        }

        @Override
        public Long call()
        {
            for(int i = 0; i < Regions.size(); ++i)
            {
                RegionData region = Regions.get(i);
                processRegion(region);

                if(i > 0 && (i % 100) == 0)
                {
                    BT_LOGGER.info("{}: processed {} regions", mTaskId, i);
                }
            }

            BT_LOGGER.debug("{}: complete", mTaskId);

            return (long)0;
        }
    }

    private void processRegion(final RegionData region)
    {
        List<RepeatMaskData> rmDataMatches = mRmAnnotation.findMatches(region);
        write(region, rmDataMatches);
    }

    private synchronized void write(final RegionData region, final List<RepeatMaskData> rmDataMatches)
    {
        try
        {
            int maxBaseOverlap = 0;
            RepeatMaskData maxOverlapRmData = null;

            for(RepeatMaskData rmData : rmDataMatches)
            {
                int baseOverlap = min(region.end(), rmData.Region.end()) - max(region.start(), rmData.Region.start()) + 1;

                if(baseOverlap > maxBaseOverlap)
                {
                    maxBaseOverlap = baseOverlap;
                    maxOverlapRmData = rmData;
                }
            }

            mWriter.write(format("%s\t%d\t%d\t%s\t%s",
                    region.Data, rmDataMatches.size(), maxBaseOverlap,
                    maxOverlapRmData != null ? maxOverlapRmData.ClassType : "", maxOverlapRmData != null ? maxOverlapRmData.Repeat : ""));

            mWriter.newLine();
         }
        catch(IOException e)
        {
            BT_LOGGER.error("failed to write annotation data: {}", e.toString());
            System.exit(1);
        }
    }

    public static void main(@NotNull final String[] args)
    {
        ConfigBuilder configBuilder = new ConfigBuilder(APP_NAME);
        configBuilder.addPath(INPUT_REGION_FILE, true, "File with enriched regions to annotation");
        configBuilder.addConfigItem(OUTPUT_FILE, false, "Output matching file");
        addRefGenomeVersion(configBuilder);
        RepeatMaskAnnotations.addConfig(configBuilder);
        addThreadOptions(configBuilder);

        addLoggingOptions(configBuilder);
        
        configBuilder.checkAndParseCommandLine(args);

        EnrichedRegionRepeatMasks enrichedRegionRepeatMasks = new EnrichedRegionRepeatMasks(configBuilder);
        enrichedRegionRepeatMasks.run();
    }
}
