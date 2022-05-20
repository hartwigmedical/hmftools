package com.hartwig.hmftools.isofox.novel.cohort;

import static com.hartwig.hmftools.common.rna.AltSpliceJunctionContext.SPLICE_JUNC;
import static com.hartwig.hmftools.common.rna.AltSpliceJunctionFile.FLD_ALT_SJ_FRAG_COUNT;
import static com.hartwig.hmftools.common.rna.AltSpliceJunctionFile.FLD_ALT_SJ_POS_END;
import static com.hartwig.hmftools.common.rna.AltSpliceJunctionFile.FLD_ALT_SJ_POS_START;
import static com.hartwig.hmftools.common.rna.AltSpliceJunctionFile.FLD_ALT_SJ_TYPE;
import static com.hartwig.hmftools.common.rna.RnaCommon.FLD_CHROMOSOME;
import static com.hartwig.hmftools.common.rna.RnaCommon.FLD_GENE_ID;
import static com.hartwig.hmftools.common.rna.RnaCommon.FLD_GENE_NAME;
import static com.hartwig.hmftools.common.utils.FileReaderUtils.createFieldsIndexMap;
import static com.hartwig.hmftools.common.utils.sv.StartEndIterator.SE_END;
import static com.hartwig.hmftools.common.utils.sv.StartEndIterator.SE_START;
import static com.hartwig.hmftools.isofox.IsofoxConfig.ISF_LOGGER;
import static com.hartwig.hmftools.isofox.cohort.AnalysisType.ALT_SPLICE_JUNCTION;
import static com.hartwig.hmftools.isofox.cohort.CohortConfig.formSampleFilenames;
import static com.hartwig.hmftools.isofox.results.ResultsWriter.DELIMITER;

import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.List;
import java.util.Map;

import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import com.hartwig.hmftools.common.rna.AltSpliceJunctionFile;
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
    private final boolean mKnownSitesOnly;

    private final AltSjWriter mWriter;

    private final AltSjFilter mAltSjFilter;

    // map of chromosomes to a map of genes to a list of alternate splice junctions
    private final Map<String,Map<String,List<AltSjCohortData>>> mAltSpliceJunctions;

    private static final String ALT_SJ_MIN_SAMPLES = "alt_sj_min_samples";
    private static final String ALT_SJ_MIN_CANCER_SAMPLES = "alt_sj_min_cancer_samples";
    private static final String ALT_SJ_PROB_THRESHOLD = "alt_sj_prob_threshold";
    private static final String ALT_SJ_MIN_FRAGS = "alt_sj_min_frags";
    private static final String ALT_SJ_KNOWN_SITES_ONLY = "alt_sj_known_sites_only";

    public AltSjCohortAnalyser(final CohortConfig config, final CommandLine cmd)
    {
        mConfig = config;
        mAltSpliceJunctions = Maps.newHashMap();

        mMinSampleThreshold = Integer.parseInt(cmd.getOptionValue(ALT_SJ_MIN_SAMPLES, "0"));
        mMinCancerSampleThreshold = Integer.parseInt(cmd.getOptionValue(ALT_SJ_MIN_CANCER_SAMPLES, "0"));
        mMinFragments = Integer.parseInt(cmd.getOptionValue(ALT_SJ_MIN_FRAGS, "0"));
        mProbabilityThreshold = Double.parseDouble(cmd.getOptionValue(ALT_SJ_PROB_THRESHOLD, "1.0"));
        mKnownSitesOnly = cmd.hasOption(ALT_SJ_KNOWN_SITES_ONLY);

        mAltSjFilter = new AltSjFilter(mConfig.RestrictedGeneIds, mConfig.ExcludedGeneIds, mMinFragments);

        boolean freqByCancerType = mConfig.SampleData.CancerTypeSamples.size() > 1 && mMinCancerSampleThreshold > 0;
        mWriter = new AltSjWriter(config, cmd, freqByCancerType);
    }

    public static void addCmdLineOptions(final Options options)
    {
        options.addOption(ALT_SJ_MIN_SAMPLES, true, "Min number of samples to report an alt SJ");
        options.addOption(ALT_SJ_MIN_CANCER_SAMPLES, true, "Min number of samples to report an alt SJ");
        options.addOption(ALT_SJ_MIN_FRAGS, true, "Min frag count supporting alt-SJs outside gene panel");
        options.addOption(ALT_SJ_PROB_THRESHOLD, true, "Only write alt SJs for fisher probability less than this");
        options.addOption(ALT_SJ_KNOWN_SITES_ONLY, false, "Only write alt SJs if at least one site is a known splice site");
        AltSjWriter.addCmdLineOptions(options);
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

                    final List<AltSpliceJunctionFile> altSJs = loadFile(altSJFile, null, mAltSjFilter);
                    ++sampleCount;

                    ISF_LOGGER.debug("{}: sample({}) loaded {} alt-SJ records", sampleCount, sampleId, altSJs.size());

                    altSJs.forEach(x -> addAltSpliceJunction(x, sampleId, cancerType));
                }

                // write out alt-SJs for this cancer type
                ISF_LOGGER.info("cancerType({}) writing alt-SJs for {} samples", cancerType, sampleIds.size());
                mWriter.writeCancerTypeFrequencies(mAltSpliceJunctions, mMinCancerSampleThreshold);
                mAltSpliceJunctions.clear();
            }
        }
        else
        {
            int totalProcessed = 0;
            int nextLog = 100000;

            // load each sample's alt SJs and consolidate into a single list
            for(int i = 0; i < mConfig.SampleData.SampleIds.size(); ++i)
            {
                final String sampleId = mConfig.SampleData.SampleIds.get(i);
                final Path altSJFile = filenames.get(i);

                final List<AltSpliceJunctionFile> altSJs = loadFile(altSJFile, null, mAltSjFilter);

                ISF_LOGGER.debug("{}: sample({}) loaded {} alt-SJ records", i, sampleId, altSJs.size());
                totalProcessed += altSJs.size();

                final String cancerType = mConfig.SampleData.SampleCancerType.get(sampleId);

                altSJs.forEach(x -> addAltSpliceJunction(x, sampleId, cancerType));

                int totalAltSJs = mAltSpliceJunctions.values().stream().mapToInt(x -> x.values().stream().mapToInt(y -> y.size()).sum()).sum();

                if(totalAltSJs >= nextLog)
                {
                    ISF_LOGGER.debug("cached alt-SJ count({})", totalAltSJs);
                    nextLog += 100000;
                }

                altSJs.forEach(x -> mWriter.writeCombinedSampleData(sampleId, x, mMinFragments));
            }

            ISF_LOGGER.info("loaded {} alt-SJ records", totalProcessed);

            if(mProbabilityThreshold < 1)
            {
                // write a report for any re-occurring alt SJ
                mWriter.writeReoccurringAltSpliceJunctions(mAltSpliceJunctions, mMinSampleThreshold, mProbabilityThreshold);
            }

            // write a cohort file
            mWriter.writeCohortFrequencies(mAltSpliceJunctions, mMinSampleThreshold);
        }

        mWriter.close();
    }

    public static List<AltSpliceJunctionFile> loadFile(final Path filename, final Map<String,Integer> refFieldsIndexMap, final AltSjFilter filter)
    {
        try
        {
            final List<String> lines = Files.readAllLines(filename);

            Map<String,Integer> fieldsIndexMap = refFieldsIndexMap != null ? refFieldsIndexMap : createFieldsIndexMap(lines.get(0), DELIMITER);

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

                try
                {
                    altSJs.add(AltSpliceJunctionFile.fromCsv(items, geneId, geneName, chr, posStart, posEnd, type,
                            fragCount, depthStart, depthEnd, regionStart, regionEnd, basesStart, basesEnd, transStart, transEnd));
                }
                catch(Exception e)
                {
                    ISF_LOGGER.error(" invalid line {}: {}", altSJs.size(), data);
                }
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
        if(mKnownSitesOnly)
        {
            if(altSJ.RegionContexts[SE_START] != SPLICE_JUNC && altSJ.RegionContexts[SE_END] != SPLICE_JUNC)
                return;
        }

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

}
