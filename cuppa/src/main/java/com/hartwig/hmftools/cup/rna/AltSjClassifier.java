package com.hartwig.hmftools.cup.rna;

import static java.lang.Math.exp;
import static java.lang.Math.log;
import static java.lang.Math.max;
import static java.lang.Math.pow;

import static com.hartwig.hmftools.common.rna.AltSpliceJunctionFile.FLD_ALT_SJ_POS_END;
import static com.hartwig.hmftools.common.rna.AltSpliceJunctionFile.FLD_ALT_SJ_POS_START;
import static com.hartwig.hmftools.common.rna.AltSpliceJunctionFile.formKey;
import static com.hartwig.hmftools.common.rna.RnaCommon.FLD_FRAG_COUNT;
import static com.hartwig.hmftools.common.stats.CosineSimilarity.calcCosineSim;
import static com.hartwig.hmftools.common.utils.file.CommonFields.FLD_CHROMOSOME;
import static com.hartwig.hmftools.common.utils.file.CommonFields.FLD_GENE_ID;
import static com.hartwig.hmftools.common.utils.file.FileDelimiters.inferFileDelimiter;
import static com.hartwig.hmftools.common.utils.file.FileWriterUtils.createBufferedReader;
import static com.hartwig.hmftools.common.utils.MatrixFile.loadMatrixDataFile;
import static com.hartwig.hmftools.common.utils.file.FileWriterUtils.closeBufferedWriter;
import static com.hartwig.hmftools.common.utils.file.FileWriterUtils.createBufferedWriter;
import static com.hartwig.hmftools.common.utils.file.FileReaderUtils.createFieldsIndexMap;
import static com.hartwig.hmftools.cup.CuppaConfig.CUP_LOGGER;
import static com.hartwig.hmftools.common.cuppa.CategoryType.ALT_SJ;
import static com.hartwig.hmftools.common.cuppa.ClassifierType.ALT_SJ_COHORT;
import static com.hartwig.hmftools.cup.common.CupConstants.ALT_SJ_DIFF_EXPONENT;
import static com.hartwig.hmftools.cup.common.CupConstants.GENE_EXP_CSS_THRESHOLD;
import static com.hartwig.hmftools.common.cuppa.ResultType.CLASSIFIER;
import static com.hartwig.hmftools.cup.common.SampleData.isKnownCancerType;
import static com.hartwig.hmftools.cup.common.SampleResult.checkIsValidCancerType;
import static com.hartwig.hmftools.cup.rna.AltSpliceJunctionPrep.loadRefAltSjIndices;
import static com.hartwig.hmftools.cup.rna.RefAltSpliceJunctions.FLD_POS_END;
import static com.hartwig.hmftools.cup.rna.RefAltSpliceJunctions.FLD_POS_START;
import static com.hartwig.hmftools.cup.rna.RefAltSpliceJunctions.loadSampleAltSjMatrixData;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.List;
import java.util.Map;

import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import com.hartwig.hmftools.common.rna.AltSpliceJunctionFile;
import com.hartwig.hmftools.common.utils.Matrix;
import com.hartwig.hmftools.common.utils.config.ConfigBuilder;
import com.hartwig.hmftools.cup.CuppaConfig;
import com.hartwig.hmftools.common.cuppa.CategoryType;
import com.hartwig.hmftools.cup.common.CuppaClassifier;
import com.hartwig.hmftools.cup.common.SampleData;
import com.hartwig.hmftools.cup.common.SampleDataCache;
import com.hartwig.hmftools.cup.common.SampleResult;
import com.hartwig.hmftools.cup.common.SampleSimilarity;

@Deprecated
public class AltSjClassifier implements CuppaClassifier
{
    private final CuppaConfig mConfig;
    private final SampleDataCache mSampleDataCache;

    // per-cancer type adjusted fragments - loads raw sample totals, calculates average per cancer type then stores log(avg + 1)
    private Matrix mRefCancerTypeMatrix;
    private final Map<String,Integer> mRefAsjIndexMap; // map from Alt-SJ into matrix rows
    private final List<String> mRefCancerTypes; // cancer types from matrix columns

    // per sample raw fragment counts per site
    private short[][] mSampleFragCounts;
    private final Map<String,Integer> mSampleIndexMap; // map from sampleId into the sample counts matrix

    private final Map<String,RnaCohortData> mCancerDataMap; // number of samples with specific read length in each cancer type

    private final double mWeightExponent;

    private BufferedWriter mCssWriter;

    private static final String LOG_CSS_VALUES = "alt_sj_log_css";
    private static final String WEIGHT_EXPONENT = "alt_sj_weight_exp";

    public AltSjClassifier(
            final CuppaConfig config, final SampleDataCache sampleDataCache, final ConfigBuilder configBuilder)
    {
        mConfig = config;
        mSampleDataCache = sampleDataCache;

        mRefAsjIndexMap = Maps.newHashMap();
        mRefCancerTypes = Lists.newArrayList();
        mRefCancerTypeMatrix = null;
        mCancerDataMap = Maps.newHashMap();

        mSampleIndexMap = Maps.newHashMap();
        mSampleFragCounts = null;
        mCssWriter = null;

        mWeightExponent = configBuilder.getDecimal(WEIGHT_EXPONENT);

        if(configBuilder.hasValue(LOG_CSS_VALUES))
        {
            initialiseCssWriter();
        }
    }

    public static void addCmdLineArgs(final ConfigBuilder configBuilder)
    {
        configBuilder.addDecimal(WEIGHT_EXPONENT, "Exponent for weighting pair-wise calcs", ALT_SJ_DIFF_EXPONENT);
        configBuilder.addFlag(LOG_CSS_VALUES, "Log CSS values");
    }

    public CategoryType categoryType() { return ALT_SJ; }
    public void close()
    {
        closeBufferedWriter(mCssWriter);
    }

    @Override
    public boolean loadData()
    {
        loadRefCancerFragCounts();

        if(mRefAsjIndexMap.isEmpty() || mRefCancerTypeMatrix == null)
            return false;

        if(mConfig.TestRefData)
        {
            if(!mConfig.RefAltSjSampleFile.isEmpty())
            {
                mSampleFragCounts = loadSampleAltSjMatrixData(mConfig.RefAltSjSampleFile, mSampleIndexMap, mRefCancerTypeMatrix.Cols);
            }
            else
            {
                CUP_LOGGER.info("missing ref cohort alt-SJ data file");
                return false;
            }
        }
        else
        {
            int sampleCount = mSampleDataCache.SampleDataList.size();
            mSampleFragCounts = new short[sampleCount][mRefCancerTypeMatrix.Cols];

            for(int i = 0; i < mSampleDataCache.SampleDataList.size(); ++i)
            {
                SampleData sample = mSampleDataCache.SampleDataList.get(i);
                loadSampleAltSJs(sample.Id, i);
                mSampleIndexMap.put(sample.Id, i);
            }
        }

        if(mSampleFragCounts == null)
        {
            return false;
        }

        if(mSampleFragCounts[0].length != mRefCancerTypeMatrix.Cols)
        {
            // keeping in mind that for the sample matrix, the sites are in the columns, whereas for the cancer ref in the rows
            CUP_LOGGER.error("alt-SJ matrix site definition mismatch: per-cancer({}) and per-sample({})",
                    mRefCancerTypeMatrix.Cols, mSampleFragCounts[0].length);
        }

        return true;
    }

    @Override
    public boolean processSample(final SampleData sample, final List<SampleResult> results, final List<SampleSimilarity> similarities)
    {
        if(mRefCancerTypeMatrix == null)
            return false;

        if(mSampleDataCache.isMultiSample() && !sample.hasRna())
            return true;

        Integer sampleIndex = mSampleIndexMap.get(sample.Id);

        if(sampleIndex == null)
        {
            CUP_LOGGER.warn("sample({}) alt-SJ frag counts matrix data not found", sample.Id);
            return false;
        }

        // prepare counts
        final short[] rawSampleFragCounts = mSampleFragCounts[sampleIndex];
        final double[] adjSampleFragCounts = convertSampleFragCounts(rawSampleFragCounts);

        addCancerCssResults(sample, rawSampleFragCounts, adjSampleFragCounts, results);

        return true;
    }

    private double[] convertSampleFragCounts(final short[] rawSampleFragCounts)
    {
        final double[] sampleFragCounts = new double[rawSampleFragCounts.length];

        for(int i = 0; i < sampleFragCounts.length; ++i)
        {
            sampleFragCounts[i] = convertFragCount(rawSampleFragCounts[i]);
        }

        return sampleFragCounts;
    }

    private void addCancerCssResults(
            final SampleData sample, final short[] rawSampleFragCounts, final double[] adjSampleFragCounts, final List<SampleResult> results)
    {
        int refCancerCount = mRefCancerTypeMatrix.Rows;

        final Map<String,Double> cancerCssTotals = Maps.newHashMap();

        int totalFrags = 0;
        int altSjSites = 0;

        if(mCssWriter != null)
        {
            for(int i = 0; i < rawSampleFragCounts.length; ++i)
            {
                totalFrags += rawSampleFragCounts[i];

                if(adjSampleFragCounts[i] > 0)
                    ++altSjSites;
            }
        }

        for(int i = 0; i < refCancerCount; ++i)
        {
            String refCancerId = mRefCancerTypes.get(i);
            final RnaCohortData cohortData = mCancerDataMap.get(refCancerId);

            if(cohortData == null)
            {
                CUP_LOGGER.error("failed to find RNA cohort data({})", refCancerId);
                continue;
            }

            if(!isKnownCancerType(cohortData.CancerType))
                continue;

            if(!checkIsValidCancerType(sample, cohortData.CancerType, cancerCssTotals))
                continue;

            if(!sample.hasRna())
                continue;

            boolean matchesCancerType = sample.cancerType().equals(cohortData.CancerType);

            final double[] refAsjFragCounts = sample.isRefSample() && matchesCancerType ?
                    adjustRefCounts(mRefCancerTypeMatrix.getRow(i), rawSampleFragCounts, cohortData.SampleCount) : mRefCancerTypeMatrix.getRow(i);

            // now any adjustments have been made, zero out any low-fragment-count sites
            double css = calcCosineSim(refAsjFragCounts, adjSampleFragCounts, false, false);

            if(css < GENE_EXP_CSS_THRESHOLD)
                continue;

            writeCssResult(sample, totalFrags, altSjSites, cohortData.CancerType, css);

            if(mSampleDataCache.isSingleSample() && CUP_LOGGER.isDebugEnabled())
            {
                int matchedSites = 0;
                int sampleSites = 0;
                int refSites = 0;

                for(int j = 0; j < adjSampleFragCounts.length; ++j)
                {
                    if(adjSampleFragCounts[j] > 0)
                        ++sampleSites;

                    if(refAsjFragCounts[j] > 0)
                        ++refSites;

                    if(adjSampleFragCounts[j] > 0 && refAsjFragCounts[j] > 0)
                        ++matchedSites;
                }

                CUP_LOGGER.debug("sample({}) cancer({}) refCancer({}) css({}) sites(sample={} ref={} matched={})",
                        sample.Id, sample.cancerType(), cohortData.CancerType, String.format("%.4f", css),
                        sampleSites, refSites, matchedSites);
            }

            double cssWeight = pow(mWeightExponent, -100 * (1 - css));

            double weightedCss = css * cssWeight;
            cancerCssTotals.put(cohortData.CancerType, weightedCss);
        }

        double totalCss = cancerCssTotals.values().stream().mapToDouble(x -> x).sum();

        for(Map.Entry<String,Double> entry : cancerCssTotals.entrySet())
        {
            cancerCssTotals.put(entry.getKey(), entry.getValue() / totalCss);
        }

        results.add(new SampleResult(
                sample.Id, ALT_SJ, CLASSIFIER, ALT_SJ_COHORT.toString(), String.format("%.4g", totalCss), cancerCssTotals));
    }

    private double[] adjustRefCounts(final double[] refCounts, final short[] sampleCounts, int cancerSampleCount)
    {
        // remove the sample's counts - which since now keep raw counts does not need to be de-logged
        double[] adjustedCounts = new double[refCounts.length];

        for(int b = 0; b < refCounts.length; ++b)
        {
            if(refCounts[b] == 0)
            {
                adjustedCounts[b] = refCounts[b];
            }
            else
            {
                double actualRefCount = (exp(refCounts[b]) - 1) * cancerSampleCount;
                double adjustedRefAvg = max(actualRefCount - sampleCounts[b], 0) / (cancerSampleCount - 1);
                adjustedCounts[b] = convertFragCount(adjustedRefAvg);
            }
        }

        return adjustedCounts;
    }

    private double convertFragCount(double fragCount)
    {
        return log(fragCount + 1);
    }

    private void loadRefCancerFragCounts()
    {
        final List<String> ignoreFields = Lists.newArrayList(FLD_GENE_ID, FLD_CHROMOSOME, FLD_POS_START, FLD_POS_END);

        if(!loadRefAltSjIndices(mConfig.RefAltSjCancerFile, mRefAsjIndexMap))
            return;

        mRefCancerTypeMatrix = loadMatrixDataFile(mConfig.RefAltSjCancerFile, mRefCancerTypes, ignoreFields, true);

        if(mRefCancerTypeMatrix == null)
            return;

        // calculate and use an average frag count per cancer type
        final double[][] refData = mRefCancerTypeMatrix.getData();
        for(int r = 0; r < mRefCancerTypeMatrix.Rows; ++r)
        {
            int cancerSampleCount = 0;
            final String cancerType = mRefCancerTypes.get(r);

            // only count those samples with RNA
            final List<SampleData> samples = mSampleDataCache.RefCancerSampleData.get(cancerType);

            if(samples != null)
            {
                cancerSampleCount = (int)samples.stream().filter(x -> x.hasRna()).count();
            }

            CUP_LOGGER.debug("alt-SJ cancerType({}) RNA samples({})", cancerType, cancerSampleCount);
            mCancerDataMap.put(cancerType, new RnaCohortData(cancerType, cancerSampleCount));

            for(int bucketIndex = 0; bucketIndex < mRefCancerTypeMatrix.Cols; ++bucketIndex)
            {
                // first calculate the average for the cancer type
                double avgFragCount = refData[r][bucketIndex] / cancerSampleCount;
                refData[r][bucketIndex] = convertFragCount(avgFragCount);
            }
        }
    }

    /*
    private boolean loadRefAltSjIndices(final String filename)
    {
        try
        {
            BufferedReader fileReader = createBufferedReader(filename);

            String header = fileReader.readLine();
            String fileDelim = inferFileDelimiter(filename);
            Map<String,Integer> fieldsIndexMap = createFieldsIndexMap(header, fileDelim);

            int chrIndex = fieldsIndexMap.get(FLD_CHROMOSOME);
            int posStartIndex = fieldsIndexMap.get(FLD_POS_START);
            int posEndIndex = fieldsIndexMap.get(FLD_POS_END);

            String line = fileReader.readLine();
            int altSjIndex = 0;

            while(line != null)
            {
                final String[] items = line.split(fileDelim, -1);

                final String asjKey = formKey(
                        items[chrIndex], Integer.parseInt(items[posStartIndex]), Integer.parseInt(items[posEndIndex]));

                mRefAsjIndexMap.put(asjKey, altSjIndex++);

                line = fileReader.readLine();
            }
        }
        catch (IOException e)
        {
            CUP_LOGGER.error("failed to read RNA ref alt-SJs from file({}): {}", filename, e.toString());
            return false;
        }

        return true;
    }
    */

    private boolean loadSampleAltSJs(final String sampleId, int sampleIndex)
    {
        final String filename = AltSpliceJunctionFile.generateFilename(mConfig.getIsofoxDataDir(sampleId), sampleId);

        if(!Files.exists(Paths.get(filename)))
            return false;

        try
        {
            final List<String> lines = Files.readAllLines(Paths.get(filename));
            String fileDelim = inferFileDelimiter(filename);

            final Map<String,Integer> fieldsIndexMap = createFieldsIndexMap(lines.get(0), fileDelim);
            lines.remove(0);

            int chrIndex = fieldsIndexMap.get(FLD_CHROMOSOME);
            int posStartIndex = fieldsIndexMap.get(FLD_ALT_SJ_POS_START);
            int posEndIndex = fieldsIndexMap.get(FLD_ALT_SJ_POS_END);
            int fragCountIndex = fieldsIndexMap.get(FLD_FRAG_COUNT);

            int matchedRefAltSJs = 0;

            for(String data : lines)
            {
                final String items[] = data.split(fileDelim, -1);

                String chromosome = items[chrIndex];
                int posStart = Integer.parseInt(items[posStartIndex]);
                int posEnd = Integer.parseInt(items[posEndIndex]);

                final String asjKey = formKey(chromosome, posStart, posEnd);

                Integer bucketIndex = mRefAsjIndexMap.get(asjKey);

                if(bucketIndex == null)
                    continue;

                short fragCount = (short)Math.min(Integer.parseInt(items[fragCountIndex]), Short.MAX_VALUE);
                ++matchedRefAltSJs;

                mSampleFragCounts[sampleIndex][bucketIndex] = fragCount;
            }

            CUP_LOGGER.info("loaded {} matching alt-SJs from file({})",
                    matchedRefAltSJs, filename);

            return matchedRefAltSJs > 0;
        }
        catch(IOException e)
        {
            CUP_LOGGER.error("failed to load alt splice junction file({}): {}", filename, e.toString());
            return false;
        }
    }

    private void initialiseCssWriter()
    {
        try
        {
            final String sampleDataFilename = mConfig.formOutputFilename("ALTSJ_CSS");

            mCssWriter = createBufferedWriter(sampleDataFilename, false);
            mCssWriter.write("SampleId,CancerType,TotalFrags,AltSjSites,RefCancerType,CSS");
            mCssWriter.newLine();
        }
        catch(IOException e)
        {
            CUP_LOGGER.error("failed to write alt-SJ CSS data: {}", e.toString());
        }
    }

    private synchronized void writeCssResult(final SampleData sample, int totalFrags, int altSjSites, final String refCancerType, double css)
    {
        if(mCssWriter == null)
            return;

        try
        {
            mCssWriter.write(String.format("%s,%s,%d,%d,%s,%.4f",
                    sample.Id, sample.cancerType(), totalFrags, altSjSites, refCancerType, css));
            mCssWriter.newLine();
        }
        catch(IOException e)
        {
            CUP_LOGGER.error("failed to write alt-SJ CSS data: {}", e.toString());
        }
    }
}
