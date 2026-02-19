package com.hartwig.hmftools.isofox.novel.cohort;

import static com.hartwig.hmftools.common.utils.file.FileDelimiters.ITEM_DELIM;
import static com.hartwig.hmftools.common.utils.file.FileWriterUtils.closeBufferedWriter;
import static com.hartwig.hmftools.common.utils.file.FileWriterUtils.createBufferedWriter;
import static com.hartwig.hmftools.common.sv.StartEndIterator.SE_END;
import static com.hartwig.hmftools.common.sv.StartEndIterator.SE_START;
import static com.hartwig.hmftools.isofox.IsofoxConfig.ISF_LOGGER;
import static com.hartwig.hmftools.isofox.novel.AltSpliceJunction.getDonorAcceptorBases;

import java.io.BufferedWriter;
import java.io.IOException;
import java.util.List;
import java.util.Map;
import java.util.StringJoiner;

import com.hartwig.hmftools.isofox.novel.AltSpliceJunctionFile;
import com.hartwig.hmftools.common.utils.config.ConfigBuilder;
import com.hartwig.hmftools.isofox.cohort.CohortConfig;

public class AltSjWriter
{
    private final CohortConfig mConfig;

    private final BufferedWriter mCombinedDataWriter; // writes alt-SJs from all samples into a single file

    private final BufferedWriter mCohortFrequencyWriter; // writes alt-SJ frequencies across all samples
    private final boolean mVerboseCohortFrequency;

    private static final String ALT_SJ_WRITE_COMBINED_COHORT = "alt_sj_write_combined_cohort";
    private static final String ALT_SJ_WRITE_COHORT_FREQ = "alt_sj_write_cohort_freq";
    private static final String ALT_SJ_COHORT_FREQ_VERBOSE = "alt_sj_cohort_freq_verbose";

    public AltSjWriter(final CohortConfig config, final ConfigBuilder configBuilder, boolean freqByCancerType)
    {
        mConfig = config;

        mCombinedDataWriter = configBuilder.hasFlag(ALT_SJ_WRITE_COMBINED_COHORT) ? initCombinedDataWriter() : null;
        mVerboseCohortFrequency = configBuilder.hasFlag(ALT_SJ_COHORT_FREQ_VERBOSE);
        mCohortFrequencyWriter = configBuilder.hasFlag(ALT_SJ_WRITE_COHORT_FREQ) ? initCohortFrequencies(freqByCancerType) : null;
    }

    public static void registerConfig(final ConfigBuilder configBuilder)
    {
        configBuilder.addFlag(ALT_SJ_WRITE_COHORT_FREQ, "Combined alt SJs from multiple samples into a single file");
        configBuilder.addFlag(ALT_SJ_WRITE_COMBINED_COHORT, "Combined alt SJs from multiple samples into a single file");
        configBuilder.addFlag(ALT_SJ_COHORT_FREQ_VERBOSE, "Combined alt SJs from multiple samples into a single file");
    }

    public void close()
    {
        closeBufferedWriter(mCombinedDataWriter);
        closeBufferedWriter(mCohortFrequencyWriter);
    }

    private BufferedWriter initCombinedDataWriter()
    {
        try
        {
            final String outputFileName = mConfig.formCohortFilename("alt_sj_cohort_sample_data.csv");
            BufferedWriter writer = createBufferedWriter(outputFileName, false);

            writer.write("SampleId,GeneId,Chromosome,Type,SjStart,SjEnd");
            writer.write(",FragCount,StartDepth,EndDepth,StartContext,EndContext,TransStart,TransEnd");
            writer.newLine();
            return writer;
        }
        catch(IOException e)
        {
            ISF_LOGGER.error("failed to initialise combined alt-SJ sample data file: {}", e.toString());
            return null;
        }
    }

    public void writeCombinedSampleData(final String sampleId, final AltSpliceJunctionFile altSJ, final int minFragments)
    {
        if(mCombinedDataWriter == null || altSJ.FragmentCount < minFragments)
            return;

        try
        {
            mCombinedDataWriter.write(String.format("%s,%s,%s,%s,%d,%d",
                    sampleId, altSJ.GeneId, altSJ.Chromosome, altSJ.Type,
                    altSJ.SpliceJunction[SE_START], altSJ.SpliceJunction[SE_END]));

            mCombinedDataWriter.write(String.format(",%d,%d,%d,%s,%s,%s,%s",
                    altSJ.FragmentCount, altSJ.DepthCounts[SE_START], altSJ.DepthCounts[SE_END],
                    altSJ.RegionContexts[SE_START], altSJ.RegionContexts[SE_END],
                    altSJ.TranscriptNames[SE_START], altSJ.TranscriptNames[SE_END]));

            mCombinedDataWriter.newLine();
        }
        catch(IOException e)
        {
            ISF_LOGGER.error("failed to write alt-SJ sample data file: {}", e.toString());
        }
    }

    private BufferedWriter initCohortFrequencies(boolean freqByCancerType)
    {
        try
        {
            final String outputFileName = mConfig.formCohortFilename("alt_sj_cohort.csv");
            BufferedWriter writer = createBufferedWriter(outputFileName, false);

            if(freqByCancerType)
            {
                writer.write("GeneId,CancerType,Chromosome,Type,SjStart,SjEnd,SampleCount,Prevalence");

                if(mVerboseCohortFrequency)
                {
                    writer.write(",SampleIds");
                }
            }
            else
            {
                writer.write("GeneId,Chromosome,Type,SjStart,SjEnd,SampleCount");

                if(mVerboseCohortFrequency)
                {
                    writer.write(",StartContext,EndContext,BaseMotif,AvgFrags,MaxFrags,AvgStartDepth,AvgEndDepth,SampleIds");
                }
            }

            writer.newLine();
            return writer;
        }
        catch(IOException e)
        {
            ISF_LOGGER.error("failed to initialise alt-SJ cohort file: {}", e.toString());
            return null;
        }
    }

    public void writeCohortFrequencies(final Map<String,Map<String, List<AltSjCohortData>>> altSpliceJunctions, final int minSampleThreshold)
    {
        if(mCohortFrequencyWriter == null)
            return;

        try
        {
            for(Map.Entry<String,Map<String, List<AltSjCohortData>>> chrEntry : altSpliceJunctions.entrySet())
            {
                final String chromosome = chrEntry.getKey();
                final Map<String,List<AltSjCohortData>> geneMap = chrEntry.getValue();

                for(Map.Entry<String,List<AltSjCohortData>> geneEntry : geneMap.entrySet())
                {
                    final String geneId = geneEntry.getKey();

                    for (AltSjCohortData altSjData : geneEntry.getValue())
                    {
                        final AltSpliceJunctionFile altSJ = altSjData.AltSJ;

                        int sampleCount = altSjData.getSampleIds().size();

                        if(sampleCount < minSampleThreshold)
                            continue;

                        mCohortFrequencyWriter.write(String.format("%s,%s,%s,%d,%d,%d",
                                geneId, chromosome, altSJ.Type, altSJ.SpliceJunction[SE_START], altSJ.SpliceJunction[SE_END], sampleCount));

                        if(mVerboseCohortFrequency)
                        {
                            mCohortFrequencyWriter.write(String.format(",%s,%s,%s",
                                    altSJ.RegionContexts[SE_START], altSJ.RegionContexts[SE_END], getDonorAcceptorBases(altSJ.BaseContexts)));

                            mCohortFrequencyWriter.write(String.format(",%.1f,%d,%.0f,%.0f",
                                    altSjData.getAvgFragmentCount(), altSjData.getMaxFragmentCount(),
                                    altSjData.getPositionCount(SE_START) / (double) sampleCount,
                                    altSjData.getPositionCount(SE_END) / (double) sampleCount));

                            mCohortFrequencyWriter.write(String.format(",%s", altSjData.sampleIdsStr()));
                        }

                        mCohortFrequencyWriter.newLine();
                    }
                }
            }
        }
        catch(IOException e)
        {
            ISF_LOGGER.error("failed to write alt-SJ cohort file: {}", e.toString());
        }
    }

    public void writeCancerTypeFrequencies(
            final Map<String,Map<String,List<AltSjCohortData>>> altSpliceJunctions, final int minCancerSampleThreshold)
    {
        if(mCohortFrequencyWriter == null)
            return;

        try
        {
            for(Map.Entry<String,Map<String,List<AltSjCohortData>>> chrEntry : altSpliceJunctions.entrySet())
            {
                final String chromosome = chrEntry.getKey();
                final Map<String,List<AltSjCohortData>> geneMap = chrEntry.getValue();

                for(Map.Entry<String,List<AltSjCohortData>> geneEntry : geneMap.entrySet())
                {
                    final String geneId = geneEntry.getKey();

                    for (AltSjCohortData altSjData : geneEntry.getValue())
                    {
                        final AltSpliceJunctionFile altSJ = altSjData.AltSJ;

                        for(Map.Entry<String,List<String>> entry : altSjData.cancerSampleIds().entrySet())
                        {
                            int sampleCount = entry.getValue().size();

                            if(sampleCount < minCancerSampleThreshold)
                                continue;

                            final String cancerType = entry.getKey();

                            int cancerSampleCount = mConfig.SampleData.CancerTypeSamples.get(cancerType).size();
                            double prevalence = sampleCount / (double)cancerSampleCount;

                            mCohortFrequencyWriter.write(String.format("%s,%s,%s,%s,%d,%d",
                                    geneId, cancerType, chromosome, altSJ.Type,
                                    altSJ.SpliceJunction[SE_START], altSJ.SpliceJunction[SE_END]));

                            mCohortFrequencyWriter.write(String.format(",%d,%.4f", sampleCount, prevalence));

                            if(mVerboseCohortFrequency)
                            {
                                StringJoiner sj = new StringJoiner(ITEM_DELIM);
                                entry.getValue().forEach(x -> sj.add(x));

                                mCohortFrequencyWriter.write(String.format(",%s", sj.toString()));
                            }

                            mCohortFrequencyWriter.newLine();
                        }
                    }
                }
            }
        }
        catch(IOException e)
        {
            ISF_LOGGER.error("failed to write alt-SJ cohort file: {}", e.toString());
        }
    }
}
