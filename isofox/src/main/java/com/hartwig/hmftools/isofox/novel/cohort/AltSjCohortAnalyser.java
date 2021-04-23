package com.hartwig.hmftools.isofox.novel.cohort;

import static com.hartwig.hmftools.common.rna.AltSpliceJunctionFile.FLD_ALT_SJ_FRAG_COUNT;
import static com.hartwig.hmftools.common.rna.AltSpliceJunctionFile.FLD_ALT_SJ_POS_END;
import static com.hartwig.hmftools.common.rna.AltSpliceJunctionFile.FLD_ALT_SJ_POS_START;
import static com.hartwig.hmftools.common.rna.AltSpliceJunctionFile.FLD_ALT_SJ_TYPE;
import static com.hartwig.hmftools.common.rna.RnaCommon.FLD_CHROMOSOME;
import static com.hartwig.hmftools.common.rna.RnaCommon.FLD_GENE_ID;
import static com.hartwig.hmftools.common.rna.RnaCommon.FLD_GENE_NAME;
import static com.hartwig.hmftools.common.utils.Strings.appendStrList;
import static com.hartwig.hmftools.common.utils.FileWriterUtils.closeBufferedWriter;
import static com.hartwig.hmftools.common.utils.FileWriterUtils.createBufferedWriter;
import static com.hartwig.hmftools.common.utils.FileWriterUtils.createFieldsIndexMap;
import static com.hartwig.hmftools.common.utils.sv.StartEndIterator.SE_END;
import static com.hartwig.hmftools.common.utils.sv.StartEndIterator.SE_START;
import static com.hartwig.hmftools.isofox.IsofoxConfig.ISF_LOGGER;
import static com.hartwig.hmftools.isofox.cohort.AnalysisType.ALT_SPLICE_JUNCTION;
import static com.hartwig.hmftools.isofox.cohort.CohortConfig.formSampleFilenames;
import static com.hartwig.hmftools.isofox.novel.AltSpliceJunction.getDonorAcceptorBases;
import static com.hartwig.hmftools.isofox.results.ResultsWriter.DELIMITER;
import static com.hartwig.hmftools.isofox.results.ResultsWriter.ITEM_DELIM;

import java.io.BufferedWriter;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.List;
import java.util.Map;

import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import com.hartwig.hmftools.common.rna.AltSpliceJunctionFile;
import com.hartwig.hmftools.common.stats.FisherExactTest;
import com.hartwig.hmftools.isofox.cohort.CohortConfig;
import com.hartwig.hmftools.isofox.cohort.SampleDataCache;

import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.Options;

public class AltSjCohortAnalyser
{
    private final CohortConfig mConfig;

    // other config
    private final int mMinSampleThreshold;
    private final int mMinCancerSampleThreshold;
    private final int mMinFragments;
    private final double mProbabilityThreshold;
    private final boolean mRewriteCohortFile;

    private BufferedWriter mCombinedDataWriter;
    private BufferedWriter mCohortDataWriter;

    private final AltSjFilter mAltSjFilter;
    private final Map<String,Integer> mFieldsMap;

    // map of chromosomes to a map of genes to a list of alternate splice junctions
    private final Map<String,Map<String,List<AltSjCohortData>>> mAltSpliceJunctions;

    private final FisherExactTest mFisherET;

    private static final String ALT_SJ_MIN_SAMPLES = "alt_sj_min_samples";
    private static final String ALT_SJ_MIN_CANCER_SAMPLES = "alt_sj_min_cancer_samples";
    private static final String ALT_SJ_PROB_THRESHOLD = "alt_sj_prob_threshold";
    private static final String ALT_SJ_MIN_FRAGS = "alt_sj_min_frags";
    private static final String ALT_SJ_REWRITE_COHORT = "alt_sj_rewrite_cohort";

    public AltSjCohortAnalyser(final CohortConfig config, final CommandLine cmd)
    {
        mConfig = config;
        mAltSpliceJunctions = Maps.newHashMap();
        mFisherET = new FisherExactTest();
        mFieldsMap = Maps.newHashMap();

        mMinSampleThreshold = Integer.parseInt(cmd.getOptionValue(ALT_SJ_MIN_SAMPLES, "0"));
        mMinCancerSampleThreshold = Integer.parseInt(cmd.getOptionValue(ALT_SJ_MIN_CANCER_SAMPLES, "0"));
        mMinFragments = Integer.parseInt(cmd.getOptionValue(ALT_SJ_MIN_FRAGS, "0"));
        mProbabilityThreshold = Double.parseDouble(cmd.getOptionValue(ALT_SJ_PROB_THRESHOLD, "1.0"));
        mRewriteCohortFile = cmd.hasOption(ALT_SJ_REWRITE_COHORT);

        mAltSjFilter = new AltSjFilter(mConfig.RestrictedGeneIds, mConfig.ExcludedGeneIds, mMinFragments);

        mCombinedDataWriter = null;
        mCohortDataWriter = null;
    }

    public static void addCmdLineOptions(final Options options)
    {
        options.addOption(ALT_SJ_MIN_SAMPLES, true, "Min number of samples to report an alt SJ");
        options.addOption(ALT_SJ_MIN_CANCER_SAMPLES, true, "Min number of samples to report an alt SJ");
        options.addOption(ALT_SJ_MIN_FRAGS, true, "Min frag count supporting alt-SJs outside gene panel");
        options.addOption(ALT_SJ_PROB_THRESHOLD, true, "Only write alt SJs for fisher probability less than this");
        options.addOption(ALT_SJ_REWRITE_COHORT, false, "Combined alt SJs from multiple samples into a single file");
    }

    public void processAltSpliceJunctions()
    {
        final List<Path> filenames = Lists.newArrayList();

        if(!formSampleFilenames(mConfig, ALT_SPLICE_JUNCTION, filenames))
            return;

        if(mConfig.SampleData.CancerTypeSamples.size() > 1 && mMinCancerSampleThreshold > 0)
        {
            int sampleCount = 0;

            // load and process samples by cancer type rather than just in order
            for(Map.Entry<String,List<String>> entry : mConfig.SampleData.CancerTypeSamples.entrySet())
            {
                final String cancerType = entry.getKey();
                final List<String> sampleIds = entry.getValue();

                ISF_LOGGER.info("cancerType({}) loading alt-SJs for {} samples", cancerType, sampleIds.size());

                for(final String sampleId : sampleIds)
                {
                    final Path altSJFile = filenames.stream().filter(x -> x.toString().contains(sampleId)).findFirst().orElse(null);

                    if(altSJFile == null)
                        continue;

                    final List<AltSpliceJunctionFile> altSJs = loadFile(altSJFile, mFieldsMap, mAltSjFilter);
                    ++sampleCount;

                    ISF_LOGGER.debug("{}: sample({}) loaded {} alt-SJ records", sampleCount, sampleId, altSJs.size());

                    altSJs.forEach(x -> addAltSpliceJunction(x, sampleId, cancerType));
                }

                // write out alt-SJs for this cancer type
                ISF_LOGGER.info("cancerType({}) writing alt-SJs for {} samples", cancerType, sampleIds.size());
                writeCancerTypeAltSpliceJunctions();
                mAltSpliceJunctions.clear();
            }
        }
        else
        {
            int totalProcessed = 0;
            int nextLog = 10000;

            // load each sample's alt SJs and consolidate into a single list
            for(int i = 0; i < mConfig.SampleData.SampleIds.size(); ++i)
            {
                final String sampleId = mConfig.SampleData.SampleIds.get(i);
                final Path altSJFile = filenames.get(i);

                final List<AltSpliceJunctionFile> altSJs = loadFile(altSJFile, mFieldsMap, mAltSjFilter);

                ISF_LOGGER.debug("{}: sample({}) loaded {} alt-SJ records", i, sampleId, altSJs.size());
                totalProcessed += altSJs.size();

                final String cancerType = mConfig.SampleData.SampleCancerType.get(sampleId);

                altSJs.forEach(x -> addAltSpliceJunction(x, sampleId, cancerType));

                int totalAltSJs = mAltSpliceJunctions.values().stream().mapToInt(x -> x.values().size()).sum();

                if(totalAltSJs >= nextLog)
                {
                    ISF_LOGGER.debug("cached alt-SJ count({})", totalAltSJs);
                    nextLog += 10000;
                }

                if(mRewriteCohortFile)
                    altSJs.forEach(x -> writeAltSpliceJunctionData(sampleId, x));
            }

            ISF_LOGGER.info("loaded {} alt-SJ records", totalProcessed);

            if(mProbabilityThreshold < 1)
            {
                // write a report for any re-occurring alt SJ
                writeReoccurringAltSpliceJunctions();
            }

            // write a cohort file
            writeCohortAltSpliceJunctions();
        }

        closeBufferedWriter(mCombinedDataWriter);
        closeBufferedWriter(mCohortDataWriter);
    }

    public static List<AltSpliceJunctionFile> loadFile(final Path filename, final Map<String,Integer> fieldsIndexMap, final AltSjFilter filter)
    {
        try
        {
            final List<String> lines = Files.readAllLines(filename);

            if(fieldsIndexMap.isEmpty())
                fieldsIndexMap.putAll(createFieldsIndexMap(lines.get(0), DELIMITER));

            lines.remove(0);

            int geneId = fieldsIndexMap.get(FLD_GENE_ID);
            int geneName = fieldsIndexMap.get(FLD_GENE_NAME);
            int chr = fieldsIndexMap.get(FLD_CHROMOSOME);
            int posStart = fieldsIndexMap.get(FLD_ALT_SJ_POS_START);
            int posEnd = fieldsIndexMap.get(FLD_ALT_SJ_POS_END);
            int type = fieldsIndexMap.get(FLD_ALT_SJ_TYPE);
            int fragCount = fieldsIndexMap.get(FLD_ALT_SJ_FRAG_COUNT);
            int depthStart = fieldsIndexMap.get("DepthStart");
            int depthEnd = fieldsIndexMap.get("DepthEnd");
            int regionStart = fieldsIndexMap.containsKey("RegionStart") ? fieldsIndexMap.get("RegionStart") : fieldsIndexMap.get("ContextStart");
            int regionEnd = fieldsIndexMap.containsKey("RegionEnd") ? fieldsIndexMap.get("RegionEnd") : fieldsIndexMap.get("ContextEnd");
            int basesStart = fieldsIndexMap.get("BasesStart");
            int basesEnd = fieldsIndexMap.get("BasesEnd");
            int transStart = fieldsIndexMap.get("TransStart");
            int transEnd = fieldsIndexMap.get("TransEnd");

            final List<AltSpliceJunctionFile> altSJs = Lists.newArrayList();

            for(String data : lines)
            {
                final String[] items = data.split(DELIMITER);

                if(!filter.passesFilter(items[geneId], Integer.parseInt(items[fragCount])))
                    continue;

                altSJs.add(AltSpliceJunctionFile.fromCsv(items, geneId, geneName, chr, posStart, posEnd, type,
                        fragCount, depthStart, depthEnd, regionStart, regionEnd, basesStart, basesEnd, transStart, transEnd));
            }

            return altSJs;
        }
        catch(IOException e)
        {
            ISF_LOGGER.error("failed to alt splice junction load file({}): {}", filename.toString(), e.toString());
            return null;
        }
    }

    private void addAltSpliceJunction(final AltSpliceJunctionFile altSJ, final String sampleId, final String cancerType)
    {
        Map<String,List<AltSjCohortData>> chrSJs = mAltSpliceJunctions.get(altSJ.Chromosome);

        if(chrSJs == null)
        {
            chrSJs = Maps.newHashMap();
            mAltSpliceJunctions.put(altSJ.Chromosome, chrSJs);
        }

        List<AltSjCohortData> geneList = chrSJs.get(altSJ.GeneId);

        if(geneList == null)
        {
            geneList = Lists.newArrayList();
            chrSJs.put(altSJ.GeneId, geneList);
        }

        AltSjCohortData asjCohortData = geneList.stream().filter(x -> x.AltSJ.matches(altSJ)).findFirst().orElse(null);

        if(asjCohortData == null)
        {
            asjCohortData = new AltSjCohortData(altSJ);
            geneList.add(asjCohortData);
        }

        if(!mConfig.SampleData.SampleCohort.isEmpty())
        {
            boolean isCohortA = mConfig.SampleData.sampleInCohort(sampleId, SampleDataCache.COHORT_A);
            asjCohortData.addSampleCount(sampleId, altSJ.FragmentCount, isCohortA);
        }
        else
        {
            asjCohortData.addSampleCount(sampleId, altSJ.FragmentCount, cancerType);
        }

        asjCohortData.addPositionCount(SE_START, altSJ.DepthCounts[SE_START]);
        asjCohortData.addPositionCount(SE_END, altSJ.DepthCounts[SE_END]);
    }

    private void writeAltSpliceJunctionData(final String sampleId, final AltSpliceJunctionFile altSJ)
    {
        if(altSJ.FragmentCount < mMinFragments)
            return;

        try
        {
            if(mCombinedDataWriter == null)
            {
                final String outputFileName = mConfig.formCohortFilename("alt_sj_cohort_sample_data.csv");
                mCombinedDataWriter = createBufferedWriter(outputFileName, false);

                mCombinedDataWriter.write("SampleId,GeneId,Chromosome,Type,SjStart,SjEnd");
                mCombinedDataWriter.write(",FragCount,StartDepth,EndDepth,StartContext,EndContext,TransStart,TransEnd");
                mCombinedDataWriter.newLine();
            }

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

    private void writeCohortAltSpliceJunctions()
    {
        try
        {
            if(mCohortDataWriter == null)
            {
                final String outputFileName = mConfig.formCohortFilename("alt_sj_cohort.csv");
                mCohortDataWriter = createBufferedWriter(outputFileName, false);

                mCohortDataWriter.write("GeneId,SampleCount,Chromosome,Type,SjStart,SjEnd");
                mCohortDataWriter.write(",StartContext,EndContext,BaseMotif,AvgFrags,MaxFrags,AvgStartDepth,AvgEndDepth,SampleIds");
                mCohortDataWriter.newLine();
            }

            for(Map.Entry<String,Map<String,List<AltSjCohortData>>> chrEntry : mAltSpliceJunctions.entrySet())
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

                        if(sampleCount < mMinSampleThreshold)
                            continue;

                        mCohortDataWriter.write(String.format("%s,%d,%s,%s,%d,%d",
                                geneId, sampleCount, chromosome, altSJ.Type, altSJ.SpliceJunction[SE_START], altSJ.SpliceJunction[SE_END]));

                        mCohortDataWriter.write(String.format(",%s,%s,%s",
                                altSJ.RegionContexts[SE_START], altSJ.RegionContexts[SE_END], getDonorAcceptorBases(altSJ.BaseContexts)));

                        mCohortDataWriter.write(String.format(",%.1f,%d,%.0f,%.0f",
                                altSjData.getAvgFragmentCount(), altSjData.getMaxFragmentCount(),
                                altSjData.getPositionCount(SE_START) / (double) sampleCount,
                                altSjData.getPositionCount(SE_END) / (double) sampleCount));

                        mCohortDataWriter.write(String.format(",%s",
                                appendStrList(altSjData.getSampleIds(), ITEM_DELIM.charAt(0))));

                        mCohortDataWriter.newLine();
                    }
                }
            }
        }
        catch(IOException e)
        {
            ISF_LOGGER.error("failed to write alt-SJ cohort file: {}", e.toString());
        }
    }

    private void writeCancerTypeAltSpliceJunctions()
    {
        try
        {
            if(mCohortDataWriter == null)
            {
                final String outputFileName = mConfig.formCohortFilename("alt_sj_cohort.csv");
                mCohortDataWriter = createBufferedWriter(outputFileName, false);

                mCohortDataWriter.write("GeneId,CancerType,Chromosome,Type,SjStart,SjEnd,SampleCount,Prevalence,SampleIds");
                mCohortDataWriter.newLine();
            }

            for(Map.Entry<String,Map<String,List<AltSjCohortData>>> chrEntry : mAltSpliceJunctions.entrySet())
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

                            if(sampleCount < mMinCancerSampleThreshold)
                                continue;

                            final String cancerType = entry.getKey();

                            int cancerSampleCount = mConfig.SampleData.CancerTypeSamples.get(cancerType).size();
                            double prevalence = sampleCount / (double)cancerSampleCount;

                            mCohortDataWriter.write(String.format("%s,%s,%s,%s,%d,%d",
                                    geneId, cancerType, chromosome, altSJ.Type,
                                    altSJ.SpliceJunction[SE_START], altSJ.SpliceJunction[SE_END]));

                            mCohortDataWriter.write(String.format(",%d,%.4f,%s",
                                    sampleCount, prevalence, appendStrList(entry.getValue(), ITEM_DELIM.charAt(0))));

                            mCohortDataWriter.newLine();
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

    private void writeReoccurringAltSpliceJunctions()
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
            mFisherET.initialise(totalSampleCount);

            int scCohortA = mConfig.SampleData.sampleCountInCohort(mConfig.SampleData.SampleIds, SampleDataCache.COHORT_A);

            for(Map.Entry<String,Map<String,List<AltSjCohortData>>> chrEntry : mAltSpliceJunctions.entrySet())
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

                        if(scWithAltSJCohortA < mMinSampleThreshold && scWithAltSJCohortB < mMinSampleThreshold)
                            continue;

                        int scNoAltSJCohortA = scCohortA - scWithAltSJCohortA;
                        int scNoAltSJCohortB = totalSampleCount - scWithAltSJCohortA - scWithAltSJCohortB - scNoAltSJCohortA;

                        double expectedVal  = scCohortA * scWithAltSJ / (double)totalSampleCount;
                        double fisherProb = mFisherET.calc(scWithAltSJCohortA, scNoAltSJCohortA, scWithAltSJCohortB, scNoAltSJCohortB, expectedVal);

                        if(fisherProb > mProbabilityThreshold)
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
