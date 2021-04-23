package com.hartwig.hmftools.cup.rna;

import static java.lang.Math.exp;
import static java.lang.Math.log;
import static java.lang.Math.max;
import static java.lang.Math.pow;
import static java.lang.Math.sqrt;

import static com.hartwig.hmftools.common.rna.AltSpliceJunctionFile.FLD_ALT_SJ_FRAG_COUNT;
import static com.hartwig.hmftools.common.rna.AltSpliceJunctionFile.FLD_ALT_SJ_POS_END;
import static com.hartwig.hmftools.common.rna.AltSpliceJunctionFile.FLD_ALT_SJ_POS_START;
import static com.hartwig.hmftools.common.rna.AltSpliceJunctionFile.formKey;
import static com.hartwig.hmftools.common.rna.RnaCommon.FLD_CHROMOSOME;
import static com.hartwig.hmftools.common.rna.RnaCommon.FLD_GENE_ID;
import static com.hartwig.hmftools.common.sigs.VectorUtils.sumVector;
import static com.hartwig.hmftools.common.stats.CosineSimilarity.calcCosineSim;
import static com.hartwig.hmftools.common.utils.MatrixUtils.loadMatrixDataFile;
import static com.hartwig.hmftools.common.utils.FileWriterUtils.closeBufferedWriter;
import static com.hartwig.hmftools.common.utils.FileWriterUtils.createBufferedWriter;
import static com.hartwig.hmftools.common.utils.FileWriterUtils.createFieldsIndexMap;
import static com.hartwig.hmftools.cup.CuppaConfig.CUP_LOGGER;
import static com.hartwig.hmftools.cup.CuppaConfig.DATA_DELIM;
import static com.hartwig.hmftools.cup.common.CategoryType.ALT_SJ;
import static com.hartwig.hmftools.cup.common.CategoryType.CLASSIFIER;
import static com.hartwig.hmftools.cup.common.ClassifierType.ALT_SJ_COHORT;
import static com.hartwig.hmftools.cup.common.ClassifierType.ALT_SJ_PAIRWISE;
import static com.hartwig.hmftools.cup.common.CupConstants.CSS_SIMILARITY_CUTOFF;
import static com.hartwig.hmftools.cup.common.CupConstants.CSS_SIMILARITY_MAX_MATCHES;
import static com.hartwig.hmftools.cup.common.CupConstants.RNA_ALT_SJ_DIFF_EXPONENT;
import static com.hartwig.hmftools.cup.common.CupConstants.RNA_GENE_EXP_CSS_THRESHOLD;
import static com.hartwig.hmftools.cup.common.ResultType.LIKELIHOOD;
import static com.hartwig.hmftools.cup.common.SampleData.RNA_READ_LENGTH_NONE;
import static com.hartwig.hmftools.cup.common.SampleData.isKnownCancerType;
import static com.hartwig.hmftools.cup.common.SampleResult.checkIsValidCancerType;
import static com.hartwig.hmftools.cup.common.SampleSimilarity.recordCssSimilarity;
import static com.hartwig.hmftools.cup.rna.RefAltSpliceJunctions.FLD_POS_END;
import static com.hartwig.hmftools.cup.rna.RefAltSpliceJunctions.FLD_POS_START;
import static com.hartwig.hmftools.cup.rna.RefAltSpliceJunctions.loadSampleAltSjMatrixData;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.FileReader;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.Arrays;
import java.util.List;
import java.util.Map;

import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import com.hartwig.hmftools.common.rna.AltSpliceJunctionFile;
import com.hartwig.hmftools.common.utils.Matrix;
import com.hartwig.hmftools.cup.CuppaConfig;
import com.hartwig.hmftools.cup.common.CategoryType;
import com.hartwig.hmftools.cup.common.CuppaClassifier;
import com.hartwig.hmftools.cup.common.SampleData;
import com.hartwig.hmftools.cup.common.SampleDataCache;
import com.hartwig.hmftools.cup.common.SampleResult;
import com.hartwig.hmftools.cup.common.SampleSimilarity;

import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.Options;

public class AltSjClassifier implements CuppaClassifier
{
    private final CuppaConfig mConfig;
    private final SampleDataCache mSampleDataCache;
    private boolean mIsValid;

    private final double mFragCountLogValue;

    // matrix CSS method
    private Matrix mRefCancerTypeMatrix;
    private final Map<String,Integer> mRefAsjIndexMap; // map from Alt-SJ into matrix rows
    private final List<String> mRefCancerTypes; // cancer types from matrix columns

    private Matrix mSampleFragCounts;
    private final Map<String,Integer> mSampleIndexMap; // map from sampleId into the sample counts matrix
    private final double mWeightExponent;
    private final boolean mRunPairwise;
    private final Map<String,RnaCohortData> mCancerDataMap; // number of samples with specific read length in each cancer type

    private BufferedWriter mCssWriter;

    private static final String SAMPLE_ALT_SJ_FILE = "sample_alt_sj_matrix_file";
    private static final String FRAG_COUNT_LOG_VALUE = "alt_sj_log_value";
    private static final String LOG_CSS_VALUES = "alt_sj_log_css";
    private static final String WEIGHT_EXPONENT = "alt_sj_weight_exp";
    private static final String RUN_PAIRWISE = "alt_sj_pairwise";

    private static final String READ_LENGTH_DELIM = "_";

    public AltSjClassifier(final CuppaConfig config, final SampleDataCache sampleDataCache, final CommandLine cmd)
    {
        mConfig = config;
        mSampleDataCache = sampleDataCache;
        mIsValid = true;

        mRefAsjIndexMap = Maps.newHashMap();
        mRefCancerTypes = Lists.newArrayList();
        mRefCancerTypeMatrix = null;
        mCancerDataMap = Maps.newHashMap();

        mSampleIndexMap = Maps.newHashMap();
        mSampleFragCounts = null;
        mCssWriter = null;

        mWeightExponent = Double.parseDouble(cmd.getOptionValue(WEIGHT_EXPONENT, String.valueOf(RNA_ALT_SJ_DIFF_EXPONENT)));

        if(cmd.hasOption(FRAG_COUNT_LOG_VALUE))
        {
            // zero will mean not applied, ie using averages as-is, otherwise value must be >= 1
            mFragCountLogValue = max(Double.parseDouble(cmd.getOptionValue(FRAG_COUNT_LOG_VALUE)), 1.0);
        }
        else
        {
            mFragCountLogValue = 1;
        }

        mRunPairwise = cmd.hasOption(RUN_PAIRWISE);

        if(config.RefAltSjCancerFile.isEmpty())
            return;

        final List<String> ignoreFields = Lists.newArrayList(FLD_GENE_ID, FLD_CHROMOSOME, FLD_POS_START, FLD_POS_END);
        loadRefFragCounts(ignoreFields);

        if(cmd.hasOption(SAMPLE_ALT_SJ_FILE))
        {
            List<String> asjLocations = Lists.newArrayList();
            mSampleFragCounts = loadSampleAltSjMatrixData(cmd.getOptionValue(SAMPLE_ALT_SJ_FILE), mSampleIndexMap, asjLocations);

            if(mSampleFragCounts ==  null)
            {
                mIsValid = false;
                return;
            }

            if(mSampleFragCounts.Cols != mRefCancerTypeMatrix.Rows)
            {
                CUP_LOGGER.error("alt-SJ matrix mismatch between ref({}) and cohort({})",
                        mRefCancerTypeMatrix.Rows, mSampleFragCounts.Rows);

                mIsValid = false;
            }

            // convert counts to log
            if(mFragCountLogValue > 0)
            {
                final double[][] sampleData = mSampleFragCounts.getData();
                for(int r = 0; r < mSampleFragCounts.Rows; ++r)
                {
                    for(int c = 0; c < mSampleFragCounts.Cols; ++c)
                    {
                        sampleData[r][c] = convertFragCount(sampleData[r][c]);
                    }
                }
            }
        }
        else
        {
            mSampleFragCounts = new Matrix(1, mRefCancerTypeMatrix.Rows);
        }

        if(cmd.hasOption(LOG_CSS_VALUES))
            initialiseCssWriter();
    }

    public static void addCmdLineArgs(Options options)
    {
        options.addOption(SAMPLE_ALT_SJ_FILE, true, "Cohort sample RNA alt-SJ frag counts matrix file");
        options.addOption(FRAG_COUNT_LOG_VALUE, true, "Use log of frag counts plus this value");
        options.addOption(WEIGHT_EXPONENT, true, "Exponent for weighting pair-wise calcs");
        options.addOption(RUN_PAIRWISE, false, "Run pair-wise classifier");
        options.addOption(LOG_CSS_VALUES, false, "Log CSS values");
    }

    public CategoryType categoryType() { return ALT_SJ; }
    public boolean isValid() { return mIsValid; }
    public void close()
    {
        closeBufferedWriter(mCssWriter);
    }

    public void processSample(final SampleData sample, final List<SampleResult> results, final List<SampleSimilarity> similarities)
    {
        if(!mIsValid || mRefCancerTypeMatrix == null)
            return;

        if(mSampleDataCache.isMultiSample() && !sample.hasRna())
            return;

        if(mSampleIndexMap.isEmpty())
        {
            if(!loadSampleAltSJsToArray(sample.Id))
                return;
        }

        Integer sampleIndex = mSampleFragCounts.Rows > 1 ? mSampleIndexMap.get(sample.Id) : 0;

        if(sampleIndex == null)
        {
            CUP_LOGGER.warn("sample({}) alt-SJ frag counts matrix data not found", sample.Id);
            return;
        }

        final double[] sampleFragCounts = mSampleFragCounts.getRow(sampleIndex);

        addCancerCssResults(sample, sampleFragCounts, results);

        if(mRunPairwise)
            addSampleCssResults(sample, sampleFragCounts, results, similarities);
    }

    private void addCancerCssResults(final SampleData sample, final double[] sampleFragCounts, final List<SampleResult> results)
    {
        int refCancerCount = mRefCancerTypeMatrix.Cols;

        final Map<String,Double> cancerCssTotals = Maps.newHashMap();

        Integer sampleIndex = mSampleFragCounts.Rows > 1 ? mSampleIndexMap.get(sample.Id) : 0;

        if(sampleIndex == null)
        {
            CUP_LOGGER.warn("sample({}) alt-SJ frag counts from matrix not found", sample.Id);
            return;
        }

        int totalFrags = (int)sumVector(sampleFragCounts);
        int altSjSites = (int)Arrays.stream(sampleFragCounts).filter(x -> x > 0).count();

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

            if(sample.rnaReadLength() != RNA_READ_LENGTH_NONE && cohortData.ReadLength != RNA_READ_LENGTH_NONE)
            {
                if(sample.rnaReadLength() != cohortData.ReadLength)
                    continue;
            }

            boolean matchesCancerType = sample.cancerType().equals(cohortData.CancerType);

            final double[] refAsjFragCounts = sample.isRefSample() && matchesCancerType ?
                    adjustRefCounts(mRefCancerTypeMatrix.getCol(i), sampleFragCounts, cohortData.SampleCount) : mRefCancerTypeMatrix.getCol(i);

            // skip any alt-SJ not present in the sample, but not the reverse
            double css = calcCosineSim(refAsjFragCounts, sampleFragCounts, false, false);

            if(css < RNA_GENE_EXP_CSS_THRESHOLD)
                continue;

            writeCssResult(sample, totalFrags, altSjSites, cohortData.CancerType, css);

            if(mSampleDataCache.isSingleSample() && CUP_LOGGER.isDebugEnabled())
            {
                int matchedSites = 0;
                int sampleSites = 0;
                int refSites = 0;

                for(int j = 0; j < mSampleFragCounts.Rows; ++j)
                {
                    if(sampleFragCounts[j] > 0)
                        ++sampleSites;

                    if(refAsjFragCounts[j] > 0)
                        ++refSites;

                    if(sampleFragCounts[j] > 0 && refAsjFragCounts[j] > 0)
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
                sample.Id, CLASSIFIER, LIKELIHOOD, ALT_SJ_COHORT.toString(), String.format("%.4g", totalCss), cancerCssTotals));
    }

    private void addSampleCssResults(
            final SampleData sample, final double[] sampleFragCounts, final List<SampleResult> results, final List<SampleSimilarity> similarities)
    {
        final Map<String,Double> cancerCssTotals = Maps.newHashMap();

        final List<SampleSimilarity> topMatches = Lists.newArrayList();

        int readLength = sample.rnaReadLength();

        for(Map.Entry<String,Integer> entry : mSampleIndexMap.entrySet())
        {
            final String refSampleId = entry.getKey();

            if(refSampleId.equals(sample.Id))
                continue;

            if(readLength != RNA_READ_LENGTH_NONE && readLength != mSampleDataCache.getRefSampleRnaReadLength(refSampleId))
                continue;

            final String refCancerType = mSampleDataCache.RefSampleCancerTypeMap.get(refSampleId);

            if(refCancerType == null)
                continue;

            if(!checkIsValidCancerType(sample, refCancerType, cancerCssTotals))
                continue;

            int refSampleCountsIndex = entry.getValue();
            final double[] refFragCount = mSampleFragCounts.getRow(refSampleCountsIndex);

            double css = calcCosineSim(sampleFragCounts, refFragCount);

            if(css < RNA_GENE_EXP_CSS_THRESHOLD)
                continue;

            if(mConfig.WriteSimilarities)
            {
                recordCssSimilarity(
                        topMatches, sample.Id, refSampleId, css, ALT_SJ_PAIRWISE.toString(),
                        CSS_SIMILARITY_MAX_MATCHES, CSS_SIMILARITY_CUTOFF);
            }

            if(!isKnownCancerType(refCancerType))
                continue;

            double cssWeight = pow(mWeightExponent, -100 * (1 - css));

            int cancerSampleCount = mSampleDataCache.getCancerSampleCount(refCancerType);
            double weightedCss = css * cssWeight / sqrt(cancerSampleCount);

            Double total = cancerCssTotals.get(refCancerType);

            if(total == null)
                cancerCssTotals.put(refCancerType, weightedCss);
            else
                cancerCssTotals.put(refCancerType, total + weightedCss);
        }

        double totalCss = cancerCssTotals.values().stream().mapToDouble(x -> x).sum();

        for(Map.Entry<String,Double> entry : cancerCssTotals.entrySet())
        {
            double prob = totalCss > 0 ? entry.getValue() / totalCss : 0;
            cancerCssTotals.put(entry.getKey(), prob);
        }

        results.add(new SampleResult(
                sample.Id, CLASSIFIER, LIKELIHOOD, ALT_SJ_PAIRWISE.toString(), String.format("%.4g", totalCss), cancerCssTotals));

        similarities.addAll(topMatches);
    }

    private double[] adjustRefCounts(final double[] refCounts, final double[] sampleCounts, int cancerSampleCount)
    {
        // remove the sample's counts after first de-logging both counts if required
        double[] adjustedCounts = new double[refCounts.length];

        for(int b = 0; b < refCounts.length; ++b)
        {
            if(mFragCountLogValue == 0)
            {
                double actualRefCount = refCounts[b] * cancerSampleCount;
                adjustedCounts[b] = max(actualRefCount - sampleCounts[b], 0) / (cancerSampleCount - 1);
            }
            else if(refCounts[b] == 0 || sampleCounts[b] == 0)
            {
                adjustedCounts[b] = refCounts[b];
            }
            else
            {
                double sampleCount = exp(sampleCounts[b]) - mFragCountLogValue;
                double actualRefCount = (exp(refCounts[b]) - mFragCountLogValue) * cancerSampleCount;
                double adjustedRefAvg = max(actualRefCount - sampleCount, 0) / (cancerSampleCount - 1);
                adjustedCounts[b] = convertFragCount(adjustedRefAvg);
            }
        }

        return adjustedCounts;
    }

    private double convertFragCount(double fragCount)
    {
        if(mFragCountLogValue < 1)
            return fragCount;

        return log(fragCount + mFragCountLogValue);
    }

    private void loadRefFragCounts(final List<String> ignoreFields)
    {
        loadRefAltSjIndices(mConfig.RefAltSjCancerFile);
        mRefCancerTypeMatrix = loadMatrixDataFile(mConfig.RefAltSjCancerFile, mRefCancerTypes, ignoreFields);

        if(mRefCancerTypeMatrix ==  null)
        {
            mIsValid = false;
            return;
        }

        // calculate and use an average frag count per cancer type
        final double[][] refData = mRefCancerTypeMatrix.getData();
        for(int c = 0; c < mRefCancerTypeMatrix.Cols; ++c)
        {
            int cancerSampleCount = 0;
            final String cancerType = mRefCancerTypes.get(c);

            if(cancerType.contains(READ_LENGTH_DELIM))
            {
                // count up the samples matching the specified read length for this cancer cohort
                final String[] cohortData = cancerType.split(READ_LENGTH_DELIM);
                final String refCancerType = cohortData[0];
                int cohortReadLength = Integer.parseInt(cohortData[1]);

                for(SampleData sample : mSampleDataCache.RefCancerSampleData.get(refCancerType))
                {
                    if(sample.rnaReadLength() != RNA_READ_LENGTH_NONE && sample.rnaReadLength() == cohortReadLength)
                        ++cancerSampleCount;
                }

                CUP_LOGGER.debug("alt-SJ cancerType({}) readLength({}) samples({})", refCancerType, cohortReadLength, cancerSampleCount);
                mCancerDataMap.put(cancerType, new RnaCohortData(refCancerType, cohortReadLength, cancerSampleCount));
            }
            else
            {
                // only count those samples with RNA
                final List<SampleData> samples = mSampleDataCache.RefCancerSampleData.get(cancerType);

                if(samples != null)
                {
                    cancerSampleCount = (int)samples.stream().filter(x -> x.rnaReadLength() != RNA_READ_LENGTH_NONE).count();
                }

                CUP_LOGGER.debug("alt-SJ cancerType({}) no specific read length, samples({})", cancerType, cancerSampleCount);
                mCancerDataMap.put(cancerType, new RnaCohortData(cancerType, RNA_READ_LENGTH_NONE, cancerSampleCount));
            }

            for(int r = 0; r < mRefCancerTypeMatrix.Rows; ++r)
            {
                // first calculate the average for the cancer type
                double avgFragCount = refData[r][c] / cancerSampleCount;
                refData[r][c] = convertFragCount(avgFragCount);
            }
        }

        mRefCancerTypeMatrix.cacheTranspose();
    }

    private boolean loadRefAltSjIndices(final String filename)
    {
        try
        {
            BufferedReader fileReader = new BufferedReader(new FileReader(filename));

            String header = fileReader.readLine();
            final Map<String,Integer> fieldsIndexMap = createFieldsIndexMap(header, DATA_DELIM);

            int chrIndex = fieldsIndexMap.get(FLD_CHROMOSOME);
            int posStartIndex = fieldsIndexMap.get(FLD_POS_START);
            int posEndIndex = fieldsIndexMap.get(FLD_POS_END);

            String line = fileReader.readLine();
            int altSjIndex = 0;

            while(line != null)
            {
                final String[] items = line.split(DATA_DELIM, -1);

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

    private boolean loadSampleAltSJsToArray(final String sampleId)
    {
        final String filename = AltSpliceJunctionFile.generateFilename(mConfig.SampleDataDir, sampleId);

        mSampleFragCounts.initialise(0);
        final double[][] sampleFragCounts = mSampleFragCounts.getData();

        if(!Files.exists(Paths.get(filename)))
            return false;

        try
        {
            final List<String> lines = Files.readAllLines(Paths.get(filename));

            final Map<String,Integer> fieldsIndexMap = createFieldsIndexMap(lines.get(0), DATA_DELIM);
            lines.remove(0);

            int chrIndex = fieldsIndexMap.get(FLD_CHROMOSOME);
            int posStartIndex = fieldsIndexMap.get(FLD_ALT_SJ_POS_START);
            int posEndIndex = fieldsIndexMap.get(FLD_ALT_SJ_POS_END);
            int fragCountIndex = fieldsIndexMap.get(FLD_ALT_SJ_FRAG_COUNT);

            int matchedRefAltSJs = 0;

            for(String data : lines)
            {
                final String items[] = data.split(DATA_DELIM, -1);

                String chromosome = items[chrIndex];
                int posStart = Integer.parseInt(items[posStartIndex]);
                int posEnd = Integer.parseInt(items[posEndIndex]);

                final String asjKey = formKey(chromosome, posStart, posEnd);

                Integer matrixIndex = mRefAsjIndexMap.get(asjKey);

                if(matrixIndex == null)
                    continue;

                int fragCount = Integer.parseInt(items[fragCountIndex]);
                ++matchedRefAltSJs;

                sampleFragCounts[0][matrixIndex] = convertFragCount(fragCount);
            }

            CUP_LOGGER.info("loaded {} matching alt-SJs from file({})", matchedRefAltSJs, filename);

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

    private void writeCssResult(final SampleData sample, int totalFrags, int altSjSites, final String refCancerType, double css)
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
