package com.hartwig.hmftools.isofox.novel.cohort;

import static com.hartwig.hmftools.common.utils.Strings.appendStrList;
import static com.hartwig.hmftools.common.utils.io.FileWriterUtils.closeBufferedWriter;
import static com.hartwig.hmftools.common.utils.io.FileWriterUtils.createBufferedWriter;
import static com.hartwig.hmftools.common.utils.io.FileWriterUtils.createFieldsIndexMap;
import static com.hartwig.hmftools.common.utils.sv.StartEndIterator.SE_END;
import static com.hartwig.hmftools.common.utils.sv.StartEndIterator.SE_START;
import static com.hartwig.hmftools.isofox.IsofoxConfig.ISF_LOGGER;
import static com.hartwig.hmftools.isofox.cohort.CohortAnalysisType.ALT_SPLICE_JUNCTION;
import static com.hartwig.hmftools.isofox.cohort.CohortConfig.formSampleFilenames;
import static com.hartwig.hmftools.isofox.novel.AltSpliceJunction.fromCsv;
import static com.hartwig.hmftools.isofox.results.ResultsWriter.DELIMITER;

import java.io.BufferedWriter;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.List;
import java.util.Map;
import java.util.stream.Collectors;

import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import com.hartwig.hmftools.common.stats.FisherExactTest;
import com.hartwig.hmftools.isofox.cohort.CohortConfig;
import com.hartwig.hmftools.isofox.cohort.SampleDataCache;
import com.hartwig.hmftools.isofox.novel.AltSpliceJunction;

import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.Options;

public class AltSpliceJunctionCohort
{
    private final CohortConfig mConfig;

    // other config
    private final int mMinSampleThreshold;
    private final int mMinFragsUnrequiredGenes;
    private final double mProbabilityThreshold;
    private final String mSpliceVariantFile;
    private final boolean mRewriteCohortFile;

    private BufferedWriter mSampleDataWriter;

    private static final String ALT_SJ_MIN_SAMPLES = "alt_sj_min_samples";
    private static final String ALT_SJ_PROB_THRESHOLD = "alt_sj_prob_threshold";
    private static final String ALT_SJ_MIN_FRAGS_REQ_GENES = "alt_sj_min_frags_req_genes";
    private static final String ALT_SJ_REWRITE_COHORT = "alt_sj_rewrite_cohort";
    private static final String SPLICE_VARIANT_FILE = "splice_variant_file";

    private final Map<String,Integer> mFieldsMap;

    // map of chromosomes to a map of genes to a list of alternate splice junctions
    private final Map<String, Map<String,List<AltSpliceJuncCohortData>>> mAltSpliceJunctions;

    private final SpliceVariantMatching mSpliceVariantMatching;

    private final FisherExactTest mFisherET;

    public AltSpliceJunctionCohort(final CohortConfig config, final CommandLine cmd)
    {
        mConfig = config;
        mAltSpliceJunctions = Maps.newHashMap();
        mFisherET = new FisherExactTest();
        mFieldsMap = Maps.newHashMap();

        mMinSampleThreshold = Integer.parseInt(cmd.getOptionValue(ALT_SJ_MIN_SAMPLES, "0"));
        mMinFragsUnrequiredGenes = Integer.parseInt(cmd.getOptionValue(ALT_SJ_MIN_FRAGS_REQ_GENES, "0"));
        mProbabilityThreshold = Double.parseDouble(cmd.getOptionValue(ALT_SJ_PROB_THRESHOLD, "1.0"));
        mSpliceVariantFile = cmd.getOptionValue(SPLICE_VARIANT_FILE);
        mRewriteCohortFile = cmd.hasOption(ALT_SJ_REWRITE_COHORT);

        mSampleDataWriter = null;

        mSpliceVariantMatching = mSpliceVariantFile != null ? new SpliceVariantMatching(mConfig, mSpliceVariantFile) : null;
    }

    public static void addCmdLineOptions(final Options options)
    {
        options.addOption(ALT_SJ_MIN_SAMPLES, true, "Min number of samples to report an alt SJ");
        options.addOption(ALT_SJ_MIN_FRAGS_REQ_GENES, true, "Min frag count supporting alt-SJs outside gene panel");
        options.addOption(ALT_SJ_PROB_THRESHOLD, true, "Only write alt SJs for fisher probability less than this");
        options.addOption(ALT_SJ_REWRITE_COHORT, false, "Combined alt SJs from multiple samples into a single file");
        options.addOption(SPLICE_VARIANT_FILE, true, "File with somatic variants potentially affecting splicing");
    }

    public void processAltSpliceJunctions()
    {
        final List<Path> filenames = Lists.newArrayList();

        if(!formSampleFilenames(mConfig, ALT_SPLICE_JUNCTION, filenames))
            return;

        int totalProcessed = 0;

        // load each sample's alt SJs and consolidate into a single list
        for(int i = 0; i < mConfig.SampleData.SampleIds.size(); ++i)
        {
            final String sampleId = mConfig.SampleData.SampleIds.get(i);
            final Path altSJFile = filenames.get(i);

            final List<AltSpliceJunction> altSJs = loadFile(altSJFile);

            ISF_LOGGER.debug("{}: sample({}) loaded {} alt-SJ records", i, sampleId, altSJs.size());
            totalProcessed += altSJs.size();

            altSJs.forEach(x -> addAltSpliceJunction(x, sampleId));

            if(mRewriteCohortFile)
                altSJs.forEach(x -> writeAltSpliceJunctionData(sampleId, x));

            if(mSpliceVariantMatching != null)
                mSpliceVariantMatching.evaluateSpliceVariants(sampleId, altSJs);
        }

        ISF_LOGGER.info("loaded {} alt-SJ records", totalProcessed);

        if(mProbabilityThreshold < 1)
        {
            // write a report for any re-occurring alt SJ
            writeReoccurringAltSpliceJunctions();
        }

        // write a cohort file
        writeCombinedAltSpliceJunctions();

        if(mSpliceVariantMatching != null)
            mSpliceVariantMatching.close();

        closeBufferedWriter(mSampleDataWriter);
    }

    private List<AltSpliceJunction> loadFile(final Path filename)
    {
        try
        {
            final List<String> lines = Files.readAllLines(filename);

            if(mFieldsMap.isEmpty())
                mFieldsMap.putAll(createFieldsIndexMap(lines.get(0), DELIMITER));

            lines.remove(0);

            return lines.stream()
                    .map(x -> fromCsv(x, mFieldsMap))
                    .filter(x -> passesFilter(x))
                    .collect(Collectors.toList());
        }
        catch(IOException e)
        {
            ISF_LOGGER.error("failed to alt splice junction load file({}): {}", filename.toString(), e.toString());
            return null;
        }
    }

    private boolean passesFilter(final AltSpliceJunction altSJ)
    {
        if(!mConfig.RestrictedGeneIds.isEmpty() && !mConfig.RestrictedGeneIds.contains(altSJ.getGeneId()))
        {
            if(mMinFragsUnrequiredGenes > 0)
            {
                if(altSJ.getFragmentCount() < mMinFragsUnrequiredGenes)
                    return false;
            }
            else
            {
                return false;
            }
        }

        if(!mConfig.ExcludedGeneIds.isEmpty() && mConfig.ExcludedGeneIds.contains(altSJ.getGeneId()))
            return false;

        return true;
    }

    private void addAltSpliceJunction(final AltSpliceJunction altSJ, final String sampleId)
    {
        Map<String,List<AltSpliceJuncCohortData>> chrSJs = mAltSpliceJunctions.get(altSJ.Chromosome);

        if(chrSJs == null)
        {
            chrSJs = Maps.newHashMap();
            mAltSpliceJunctions.put(altSJ.Chromosome, chrSJs);
        }

        List<AltSpliceJuncCohortData> geneList = chrSJs.get(altSJ.getGeneId());

        if(geneList == null)
        {
            geneList = Lists.newArrayList();
            chrSJs.put(altSJ.getGeneId(), geneList);
        }

        AltSpliceJuncCohortData altSjData = geneList.stream().filter(x -> x.AltSJ.matches(altSJ)).findFirst().orElse(null);

        if(altSjData == null)
        {
            altSjData = new AltSpliceJuncCohortData(altSJ);
            geneList.add(altSjData);
        }

        if(!mConfig.SampleData.SampleCohort.isEmpty())
        {
            boolean isCohortA = mConfig.SampleData.sampleInCohort(sampleId, SampleDataCache.COHORT_A);
            altSjData.addSampleAndCount(sampleId, altSJ.getFragmentCount(), isCohortA);
        }
        else
        {
            altSjData.addSampleAndCount(sampleId, altSJ.getFragmentCount());
        }

        altSjData.addPositionCount(SE_START, altSJ.getPositionCount(SE_START));
        altSjData.addPositionCount(SE_END, altSJ.getPositionCount(SE_END));
    }

    private void writeAltSpliceJunctionData(final String sampleId, final AltSpliceJunction altSJ)
    {
        if(altSJ.getFragmentCount() < mMinFragsUnrequiredGenes)
            return;

        try
        {
            if(mSampleDataWriter == null)
            {
                final String outputFileName = mConfig.formCohortFilename("alt_sj_cohort_sample_data.csv");
                mSampleDataWriter = createBufferedWriter(outputFileName, false);

                mSampleDataWriter.write("SampleId,GeneId,Chromosome,Type,SjStart,SjEnd");
                mSampleDataWriter.write(",FragCount,StartDepth,EndDepth,StartContext,EndContext,TransStart,TransEnd");
                mSampleDataWriter.newLine();
            }

            mSampleDataWriter.write(String.format("%s,%s,%s,%s,%d,%d",
                    sampleId, altSJ.getGeneId(), altSJ.Chromosome, altSJ.type(),
                    altSJ.SpliceJunction[SE_START], altSJ.SpliceJunction[SE_END]));

            mSampleDataWriter.write(String.format(",%d,%d,%d,%s,%s,%s,%s",
                    altSJ.getFragmentCount(), altSJ.getPositionCount(SE_START), altSJ.getPositionCount(SE_END),
                    altSJ.RegionContexts[SE_START], altSJ.RegionContexts[SE_END],
                    altSJ.getTranscriptNames()[SE_START], altSJ.getTranscriptNames()[SE_END]));

            mSampleDataWriter.newLine();
        }
        catch(IOException e)
        {
            ISF_LOGGER.error("failed to write alt-SJ sample data file: {}", e.toString());
        }
    }

    private void writeCombinedAltSpliceJunctions()
    {
        try
        {
            final String outputFileName = mConfig.formCohortFilename("alt_sj_cohort.csv");
            final BufferedWriter writer = createBufferedWriter(outputFileName, false);

            writer.write("GeneId,Chromosome,Type,SjStart,SjEnd,SampleCount");
            writer.write(",StartContext,EndContext,BaseMotif");
            writer.write(",AvgFrags,MaxFrags,AvgStartDepth,AvgEndDepth,SampleIds");
            writer.newLine();

            for(Map.Entry<String,Map<String,List<AltSpliceJuncCohortData>>> chrEntry : mAltSpliceJunctions.entrySet())
            {
                final String chromosome = chrEntry.getKey();
                final Map<String,List<AltSpliceJuncCohortData>> geneMap = chrEntry.getValue();

                for(Map.Entry<String,List<AltSpliceJuncCohortData>> geneEntry : geneMap.entrySet())
                {
                    final String geneId = geneEntry.getKey();

                    for (AltSpliceJuncCohortData altSjData : geneEntry.getValue())
                    {
                        final AltSpliceJunction altSJ = altSjData.AltSJ;

                        int sampleCount = altSjData.getSampleIds().size();

                        if(sampleCount < mMinSampleThreshold)
                            continue;

                        writer.write(String.format("%s,%s,%s,%d,%d,%d",
                                geneId, chromosome, altSJ.type(), altSJ.SpliceJunction[SE_START], altSJ.SpliceJunction[SE_END], sampleCount));

                        writer.write(String.format(",%s,%s,%s",
                                altSJ.RegionContexts[SE_START], altSJ.RegionContexts[SE_END], altSJ.getDonorAcceptorBases()));

                        writer.write(String.format(",%.0f,%d,%.0f,%.0f",
                                altSjData.getAvgFragmentCount(), altSjData.getMaxFragmentCount(),
                                altSjData.getPositionCount(SE_START)/(double)sampleCount,
                                altSjData.getPositionCount(SE_END)/(double)sampleCount));

                        writer.write(String.format(",%s",
                                appendStrList(altSjData.getSampleIds(), ';')));

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

            for(Map.Entry<String,Map<String,List<AltSpliceJuncCohortData>>> chrEntry : mAltSpliceJunctions.entrySet())
            {
                final String chromosome = chrEntry.getKey();
                final Map<String,List<AltSpliceJuncCohortData>> geneMap = chrEntry.getValue();

                for(Map.Entry<String,List<AltSpliceJuncCohortData>> geneEntry : geneMap.entrySet())
                {
                    final String geneId = geneEntry.getKey();

                    for (AltSpliceJuncCohortData altSjData : geneEntry.getValue())
                    {
                        final AltSpliceJunction altSJ = altSjData.AltSJ;

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
                                geneId, chromosome, altSJ.type(),
                                altSJ.SpliceJunction[SE_START], altSJ.SpliceJunction[SE_END]));

                        writer.write(String.format(",%s,%s,%s",
                                altSJ.RegionContexts[SE_START], altSJ.RegionContexts[SE_END], altSJ.getDonorAcceptorBases()));

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
