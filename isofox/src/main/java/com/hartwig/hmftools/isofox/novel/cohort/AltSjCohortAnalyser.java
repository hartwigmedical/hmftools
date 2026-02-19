package com.hartwig.hmftools.isofox.novel.cohort;

import static com.hartwig.hmftools.common.rna.AltSpliceJunctionContext.SPLICE_JUNC;
import static com.hartwig.hmftools.common.rna.NovelSpliceJunctionFile.FLD_ALT_SJ_POS_END;
import static com.hartwig.hmftools.common.rna.NovelSpliceJunctionFile.FLD_ALT_SJ_POS_START;
import static com.hartwig.hmftools.common.rna.NovelSpliceJunctionFile.FLD_ALT_SJ_TYPE;
import static com.hartwig.hmftools.common.rna.NovelSpliceJunctionFile.FLD_BASES_END;
import static com.hartwig.hmftools.common.rna.NovelSpliceJunctionFile.FLD_BASES_START;
import static com.hartwig.hmftools.common.rna.NovelSpliceJunctionFile.FLD_COHORT_FREQUENCY;
import static com.hartwig.hmftools.common.utils.file.CommonFields.FLD_REGION_END;
import static com.hartwig.hmftools.common.utils.file.CommonFields.FLD_REGION_START;
import static com.hartwig.hmftools.common.utils.file.FileDelimiters.TSV_DELIM;
import static com.hartwig.hmftools.common.utils.file.CommonFields.FLD_CHROMOSOME;
import static com.hartwig.hmftools.common.rna.RnaCommon.FLD_DEPTH_END;
import static com.hartwig.hmftools.common.rna.RnaCommon.FLD_DEPTH_START;
import static com.hartwig.hmftools.common.rna.RnaCommon.FLD_FRAG_COUNT;
import static com.hartwig.hmftools.common.utils.file.CommonFields.FLD_GENE_ID;
import static com.hartwig.hmftools.common.utils.file.CommonFields.FLD_GENE_NAME;
import static com.hartwig.hmftools.common.utils.file.FileDelimiters.inferFileDelimiter;
import static com.hartwig.hmftools.common.utils.file.FileReaderUtils.createFieldsIndexMap;
import static com.hartwig.hmftools.common.sv.StartEndIterator.SE_END;
import static com.hartwig.hmftools.common.sv.StartEndIterator.SE_START;
import static com.hartwig.hmftools.isofox.IsofoxConfig.ISF_LOGGER;
import static com.hartwig.hmftools.isofox.cohort.AnalysisType.ALT_SPLICE_JUNCTION;
import static com.hartwig.hmftools.isofox.cohort.AnalysisType.CANONICAL_SPLICE_JUNCTION;
import static com.hartwig.hmftools.isofox.cohort.CohortConfig.formSampleFilenames;
import static com.hartwig.hmftools.isofox.novel.AltSpliceJunctionFile.FLD_TRANS_END;
import static com.hartwig.hmftools.isofox.novel.AltSpliceJunctionFile.FLD_TRANS_START;

import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.List;
import java.util.Map;

import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import com.hartwig.hmftools.common.rna.AltSpliceJunctionContext;
import com.hartwig.hmftools.isofox.novel.AltSpliceJunctionFile;
import com.hartwig.hmftools.common.rna.AltSpliceJunctionType;
import com.hartwig.hmftools.common.utils.config.ConfigBuilder;
import com.hartwig.hmftools.isofox.cohort.AnalysisType;
import com.hartwig.hmftools.isofox.cohort.CohortConfig;

public class AltSjCohortAnalyser
{
    private final CohortConfig mConfig;

    // other config
    private final int mMinSampleThreshold;
    private final int mMinCancerSampleThreshold;
    private final int mMinFragments;
    private final boolean mKnownSitesOnly;
    private final boolean mLoadCanonical;

    private final AltSjWriter mWriter;

    private final AltSjFilter mAltSjFilter;

    // map of chromosomes to a map of genes to a list of alternate splice junctions
    private final Map<String,Map<String,List<AltSjCohortData>>> mAltSpliceJunctions;

    private static final String ALT_SJ_MIN_SAMPLES = "alt_sj_min_samples";
    private static final String ALT_SJ_MIN_CANCER_SAMPLES = "alt_sj_min_cancer_samples";
    private static final String ALT_SJ_MIN_FRAGS = "alt_sj_min_frags";
    private static final String ALT_SJ_KNOWN_SITES_ONLY = "alt_sj_known_sites_only";
    public static final String ALT_SJ_LOAD_CANONICAL = "alt_sj_load_canonical";
    public static final String ALT_SJ_LOAD_CANONICAL_DESC = "Load canonical splice junction files isntead of alts";

    public AltSjCohortAnalyser(final CohortConfig config, final ConfigBuilder configBuilder)
    {
        mConfig = config;
        mAltSpliceJunctions = Maps.newHashMap();

        mMinSampleThreshold = configBuilder.getInteger(ALT_SJ_MIN_SAMPLES);
        mMinCancerSampleThreshold = configBuilder.getInteger(ALT_SJ_MIN_CANCER_SAMPLES);
        mMinFragments = configBuilder.getInteger(ALT_SJ_MIN_FRAGS);
        mLoadCanonical = configBuilder.hasFlag(ALT_SJ_LOAD_CANONICAL);
        mKnownSitesOnly = configBuilder.hasFlag(ALT_SJ_KNOWN_SITES_ONLY);

        mAltSjFilter = new AltSjFilter(mConfig.RestrictedGeneIds, mConfig.ExcludedGeneIds, mMinFragments);

        boolean freqByCancerType = mConfig.SampleData.CancerTypeSamples.size() > 1 && mMinCancerSampleThreshold > 0;
        mWriter = new AltSjWriter(config, configBuilder, freqByCancerType);
    }

    public static void registerConfig(final ConfigBuilder configBuilder)
    {
        configBuilder.addInteger(ALT_SJ_MIN_SAMPLES, "Min number of samples to report an alt SJ", 0);
        configBuilder.addInteger(ALT_SJ_MIN_CANCER_SAMPLES, "Min number of samples to report an alt SJ", 0);
        configBuilder.addInteger(ALT_SJ_MIN_FRAGS, "Min frag count supporting alt-SJs outside gene panel", 0);
        configBuilder.addFlag(ALT_SJ_LOAD_CANONICAL, ALT_SJ_LOAD_CANONICAL_DESC);
        configBuilder.addFlag(ALT_SJ_KNOWN_SITES_ONLY, "Only write alt SJs if at least one site is a known splice site");
        AltSjWriter.registerConfig(configBuilder);
    }

    public void processAltSpliceJunctions()
    {
        final List<Path> filenames = Lists.newArrayList();

        AnalysisType analysisType = mLoadCanonical ? CANONICAL_SPLICE_JUNCTION : ALT_SPLICE_JUNCTION;
        if(!formSampleFilenames(mConfig, analysisType, filenames))
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

                    final List<AltSpliceJunctionFile> altSJs = mLoadCanonical ?
                            loadCanonicalSpliceFile(altSJFile, null, mAltSjFilter) :
                            loadFile(altSJFile, null, mAltSjFilter);

                    ++sampleCount;

                    ISF_LOGGER.debug("{}: sample({}) loaded {} splice-junction records", sampleCount, sampleId, altSJs.size());

                    altSJs.forEach(x -> addAltSpliceJunction(x, sampleId, cancerType));
                }

                // write out alt-SJs for this cancer type
                ISF_LOGGER.info("cancerType({}) writing splice-junction for {} samples", cancerType, sampleIds.size());
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

            // write a cohort file
            mWriter.writeCohortFrequencies(mAltSpliceJunctions, mMinSampleThreshold);
        }

        mWriter.close();
    }

    public static List<AltSpliceJunctionFile> loadFile(final Path filename, final Map<String,Integer> refFieldsIndexMap, final AltSjFilter filter)
    {
        try
        {
            List<String> lines = Files.readAllLines(filename);
            String fileDelim = inferFileDelimiter(filename.toString());

            Map<String,Integer> fieldsIndexMap = refFieldsIndexMap != null ? refFieldsIndexMap : createFieldsIndexMap(lines.get(0), fileDelim);

            lines.remove(0);

            int geneIdIndex = fieldsIndexMap.get(FLD_GENE_ID);
            int geneNameIndex = fieldsIndexMap.get(FLD_GENE_NAME);
            int chrIndex = fieldsIndexMap.get(FLD_CHROMOSOME);
            int posStartIndex = fieldsIndexMap.get(FLD_ALT_SJ_POS_START);
            int posEndIndex = fieldsIndexMap.get(FLD_ALT_SJ_POS_END);
            int typeIndex = fieldsIndexMap.get(FLD_ALT_SJ_TYPE);
            int fragCountIndex = fieldsIndexMap.get(FLD_FRAG_COUNT);
            int depthStartIndex = fieldsIndexMap.get(FLD_DEPTH_START);
            int depthEndIndex = fieldsIndexMap.get(FLD_DEPTH_END);
            int regionStartIndex = fieldsIndexMap.get(FLD_REGION_START);
            int regionEndIndex = fieldsIndexMap.get(FLD_REGION_END);
            int basesStartIndex = fieldsIndexMap.get(FLD_BASES_START);
            int basesEndIndex = fieldsIndexMap.get(FLD_BASES_END);
            int transStartIndex = fieldsIndexMap.get(FLD_TRANS_START);
            int transEndIndex = fieldsIndexMap.get(FLD_TRANS_END);
            int cohortFreqIndex = fieldsIndexMap.get(FLD_COHORT_FREQUENCY);

            final List<AltSpliceJunctionFile> altSJs = Lists.newArrayList();

            for(String data : lines)
            {
                String[] values = data.split(fileDelim);

                if(!filter.passesFilter(values[geneIdIndex], Integer.parseInt(values[fragCountIndex])))
                    continue;

                try
                {
                    altSJs.add(parseLine(values, geneIdIndex, geneNameIndex, chrIndex, posStartIndex, posEndIndex, typeIndex,
                            fragCountIndex, depthStartIndex, depthEndIndex, regionStartIndex, regionEndIndex, basesStartIndex, basesEndIndex,
                            transStartIndex, transEndIndex, cohortFreqIndex));
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
            ISF_LOGGER.error("failed to load alt splice junction file({}): {}", filename.toString(), e.toString());
            return null;
        }
    }

    private static AltSpliceJunctionFile parseLine(
            final String[] values,
            int geneId, int geneName, int chr, int posStart, int posEnd, int type, int fragCount, int depthStart, int depthEnd,
            int regionStart, int regionEnd, int basesStart, int basesEnd, int transStart, int transEnd, int cohortFrequency)
    {
        return new AltSpliceJunctionFile(
                values[geneId], values[geneName], values[chr], new int[] { Integer.parseInt(values[posStart]), Integer.parseInt(values[posEnd]) },
                AltSpliceJunctionType.valueOf(values[type]),
                Integer.parseInt(values[fragCount]), new int[] { Integer.parseInt(values[depthStart]), Integer.parseInt(values[depthEnd]) },
                new AltSpliceJunctionContext[] { AltSpliceJunctionContext.valueOf(values[regionStart]), AltSpliceJunctionContext.valueOf(values[regionEnd]) },
                new String[] { values[basesStart], values[basesEnd] }, new String[] { values[transStart], values[transEnd]},
                Integer.parseInt(values[cohortFrequency]));
    }


    public static List<AltSpliceJunctionFile> loadCanonicalSpliceFile(
            final Path filename, final Map<String,Integer> refFieldsIndexMap, final AltSjFilter filter)
    {
        try
        {
            List<String> lines = Files.readAllLines(filename);

            String fileDelim = inferFileDelimiter(filename.toString());
            Map<String,Integer> fieldsIndexMap = refFieldsIndexMap != null ? refFieldsIndexMap : createFieldsIndexMap(lines.get(0), fileDelim);

            lines.remove(0);

            int geneIdIndex = fieldsIndexMap.get(FLD_GENE_ID);
            int geneNameIndex = fieldsIndexMap.get(FLD_GENE_NAME);
            int chrIndex = fieldsIndexMap.get(FLD_CHROMOSOME);
            int posStartIndex = fieldsIndexMap.get(FLD_ALT_SJ_POS_START);
            int posEndIndex = fieldsIndexMap.get(FLD_ALT_SJ_POS_END);
            int fragCountIndex = fieldsIndexMap.get(FLD_FRAG_COUNT);
            int depthStartIndex = fieldsIndexMap.get(FLD_DEPTH_START);
            int depthEndIndex = fieldsIndexMap.get(FLD_DEPTH_END);

            final List<AltSpliceJunctionFile> altSJs = Lists.newArrayList();

            final AltSpliceJunctionContext[] emptyRegionContexts = { SPLICE_JUNC, SPLICE_JUNC };
            final String[] emptyBaseContexts = {"", ""};
            final String[] emptyTranscriptNames = {"", ""};

            for(String data : lines)
            {
                String[] values = data.split(fileDelim);

                String geneId = values[geneIdIndex];
                int fragCount = Integer.parseInt(values[fragCountIndex]);

                if(!filter.passesFilter(geneId, fragCount))
                    continue;

                int[] spliceJunction = new int[] { Integer.parseInt(values[posStartIndex]), Integer.parseInt(values[posEndIndex]) };
                int[] depthCounts = new int[] { Integer.parseInt(values[depthStartIndex]), Integer.parseInt(values[depthEndIndex]) };

                altSJs.add(new AltSpliceJunctionFile(
                        geneId, values[geneNameIndex], values[chrIndex], spliceJunction, AltSpliceJunctionType.CANONICAL,
                        fragCount, depthCounts, emptyRegionContexts, emptyBaseContexts, emptyTranscriptNames, 0));
            }

            return altSJs;
        }
        catch(IOException e)
        {
            ISF_LOGGER.error("failed to load canonical splice junction file({}): {}", filename.toString(), e.toString());
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

        asjCohortData.addSampleCount(sampleId, altSJ.FragmentCount, cancerType);

        asjCohortData.addPositionCount(SE_START, altSJ.DepthCounts[SE_START]);
        asjCohortData.addPositionCount(SE_END, altSJ.DepthCounts[SE_END]);
    }

}
