package com.hartwig.hmftools.isofox.novel.cohort;

import static com.hartwig.hmftools.common.utils.FileWriterUtils.closeBufferedWriter;
import static com.hartwig.hmftools.common.utils.FileWriterUtils.createBufferedWriter;
import static com.hartwig.hmftools.common.utils.Strings.appendStrList;
import static com.hartwig.hmftools.common.utils.sv.StartEndIterator.SE_END;
import static com.hartwig.hmftools.common.utils.sv.StartEndIterator.SE_START;
import static com.hartwig.hmftools.isofox.IsofoxConfig.ISF_LOGGER;
import static com.hartwig.hmftools.isofox.novel.AltSpliceJunction.getDonorAcceptorBases;
import static com.hartwig.hmftools.isofox.results.ResultsWriter.ITEM_DELIM;

import java.io.BufferedWriter;
import java.io.IOException;
import java.util.List;
import java.util.Map;

import com.hartwig.hmftools.common.rna.AltSpliceJunctionFile;
import com.hartwig.hmftools.common.stats.FisherExactTest;
import com.hartwig.hmftools.isofox.cohort.CohortConfig;
import com.hartwig.hmftools.isofox.cohort.SampleDataCache;

import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.Options;

public class AltSjWriter
{
    private final CohortConfig mConfig;

    private final BufferedWriter mCombinedDataWriter; // writes alt-SJs from all samples into a single file

    private final BufferedWriter mCohortFrequencyWriter; // writes alt-SJ frequencies across all samples
    private final boolean mVerboseCohortFrequency;

    private static final String ALT_SJ_WRITE_COMBINED_COHORT = "alt_sj_write_combined_cohort";
    private static final String ALT_SJ_WRITE_COHORT_FREQ = "alt_sj_write_cohort_freq";
    private static final String ALT_SJ_COHORT_FREQ_VERBOSE = "alt_sj_cohort_freq_verbose";

    public AltSjWriter(final CohortConfig config, final CommandLine cmd, boolean freqByCancerType)
    {
        mConfig = config;

        mCombinedDataWriter = cmd.hasOption(ALT_SJ_WRITE_COMBINED_COHORT) ? initCombinedDataWriter() : null;

        mVerboseCohortFrequency = cmd.hasOption(ALT_SJ_COHORT_FREQ_VERBOSE);
        mCohortFrequencyWriter = cmd.hasOption(ALT_SJ_WRITE_COHORT_FREQ) ? initCohortFrequencies(freqByCancerType) : null;
    }

    public static void addCmdLineOptions(final Options options)
    {
        options.addOption(ALT_SJ_WRITE_COHORT_FREQ, false, "Combined alt SJs from multiple samples into a single file");
        options.addOption(ALT_SJ_WRITE_COMBINED_COHORT, false, "Combined alt SJs from multiple samples into a single file");
        options.addOption(ALT_SJ_COHORT_FREQ_VERBOSE, false, "Combined alt SJs from multiple samples into a single file");
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
                writer.write("GeneId,CancerType,Chromosome,Type,SjStart,SjEnd,SampleCount,Prevalence,SampleIds");
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

                            mCohortFrequencyWriter.write(String.format(",%d,%.4f,%s",
                                    sampleCount, prevalence, appendStrList(entry.getValue(), ITEM_DELIM.charAt(0))));

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

    public void writeReoccurringAltSpliceJunctions(
            final Map<String,Map<String,List<AltSjCohortData>>> altSpliceJunctions,
            final int minSampleThreshold, final double probabilityThreshold)
    {
        try
        {
            final String outputFileName = mConfig.formCohortFilename("alt_sj_cohort_compare.csv");
            final BufferedWriter writer = createBufferedWriter(outputFileName, false);

            writer.write("GeneId,Chromosome,Type,SjStart,SjEnd");
            writer.write(",StartContext,EndContext,BaseMotif,AvgFragsCohortA,AvgFragsCohortB,MaxFragsCohortA,MaxFragsCohortB");
            writer.write(",AvgStartDepth,AvgEndDepth,AltSJSampleCount,FetProb,ExpVal");
            writer.write(",CohortAWithAltSJ,CohortBWithAltSJ,CohortA,CohortANoAltSJ,CohortBNoAltSJ");
            writer.write(",CohortASampleIds,CohortBSampleIds");
            writer.newLine();

            int totalSampleCount = mConfig.SampleData.SampleIds.size();

            final FisherExactTest fisherET = new FisherExactTest();
            fisherET.initialise(totalSampleCount);

            int scCohortA = mConfig.SampleData.sampleCountInCohort(mConfig.SampleData.SampleIds, SampleDataCache.COHORT_A);

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

                        int scWithAltSJCohortA = altSjData.getCohortSampleIds(true).size();
                        int scWithAltSJCohortB = altSjData.getCohortSampleIds(false).size();
                        int scWithAltSJ = scWithAltSJCohortA + scWithAltSJCohortB;

                        if(scWithAltSJCohortA < minSampleThreshold && scWithAltSJCohortB < minSampleThreshold)
                            continue;

                        int scNoAltSJCohortA = scCohortA - scWithAltSJCohortA;
                        int scNoAltSJCohortB = totalSampleCount - scWithAltSJCohortA - scWithAltSJCohortB - scNoAltSJCohortA;

                        double expectedVal  = scCohortA * scWithAltSJ / (double)totalSampleCount;
                        double fisherProb = fisherET.calc(scWithAltSJCohortA, scNoAltSJCohortA, scWithAltSJCohortB, scNoAltSJCohortB, expectedVal);

                        if(fisherProb > probabilityThreshold)
                            continue;

                        writer.write(String.format("%s,%s,%s,%d,%d",
                                geneId, chromosome, altSJ.Type,
                                altSJ.SpliceJunction[SE_START], altSJ.SpliceJunction[SE_END]));

                        writer.write(String.format(",%s,%s,%s",
                                altSJ.RegionContexts[SE_START], altSJ.RegionContexts[SE_END], getDonorAcceptorBases(altSJ.BaseContexts)));

                        writer.write(String.format(",%.0f,%.0f,%d,%d,%.0f,%.0f",
                                altSjData.getAvgFragmentCount(true), altSjData.getAvgFragmentCount(false),
                                altSjData.getMaxFragmentCount(true), altSjData.getMaxFragmentCount(false),
                                altSjData.getPositionCount(SE_START)/(double)scWithAltSJ,
                                altSjData.getPositionCount(SE_END)/(double)scWithAltSJ));

                        writer.write(String.format(",%d,%4.3e,%.1f,%d,%d,%d,%d,%d",
                                scWithAltSJ, fisherProb, expectedVal,
                                scWithAltSJCohortA, scWithAltSJCohortB, scCohortA, scNoAltSJCohortA, scNoAltSJCohortB));

                        writer.write(String.format(",%s,%s",
                                appendStrList(altSjData.getCohortSampleIds(true), ';'),
                                appendStrList(altSjData.getCohortSampleIds(false), ';')));

                        writer.newLine();
                    }
                }
            }

            closeBufferedWriter(writer);
        }
        catch(IOException e)
        {
            ISF_LOGGER.error("failed to write alt-SJ cohort file: {}", e.toString());
        }
    }

}
