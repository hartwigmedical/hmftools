package com.hartwig.hmftools.cup.rna;

import static java.lang.Math.exp;
import static java.lang.Math.log;
import static java.lang.Math.max;
import static java.lang.Math.pow;

import static com.hartwig.hmftools.common.stats.CosineSimilarity.calcCosineSim;
import static com.hartwig.hmftools.common.utils.MatrixUtils.loadMatrixDataFile;
import static com.hartwig.hmftools.common.utils.io.FileWriterUtils.createFieldsIndexMap;
import static com.hartwig.hmftools.cup.CuppaConfig.CUP_LOGGER;
import static com.hartwig.hmftools.cup.CuppaConfig.DATA_DELIM;
import static com.hartwig.hmftools.cup.common.CategoryType.ALT_SJ;
import static com.hartwig.hmftools.cup.common.CategoryType.CLASSIFIER;
import static com.hartwig.hmftools.cup.common.CupCalcs.adjustLowProbabilities;
import static com.hartwig.hmftools.cup.common.CupConstants.RNA_ALT_SJ_DIFF_EXPONENT;
import static com.hartwig.hmftools.cup.common.CupConstants.RNA_GENE_EXP_CSS_THRESHOLD;
import static com.hartwig.hmftools.cup.common.ResultType.LIKELIHOOD;
import static com.hartwig.hmftools.cup.common.SampleData.isKnownCancerType;
import static com.hartwig.hmftools.cup.common.SampleResult.checkIsValidCancerType;

import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.List;
import java.util.Map;
import java.util.Set;

import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import com.google.common.collect.Sets;
import com.hartwig.hmftools.common.utils.Matrix;
import com.hartwig.hmftools.cup.CuppaConfig;
import com.hartwig.hmftools.cup.common.CategoryType;
import com.hartwig.hmftools.cup.common.ClassifierType;
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

    // prevalence method
    private final Map<String,List<AltSjPrevData>> mRefAltSjPrevalence; // ref alt-SJs by chromosome

    // matrix CSS method
    private Matrix mRefCancerTypeMatrix;
    private final Map<String,Integer> mRefAsjIndexMap; // map from Alt-SJ into matrix rows
    private final List<String> mRefCancerTypes; // cancer types from matrix columns
    private final Set<Integer> mRefAltSjStartPositions; // for fast sample data filtering

    private Matrix mSampleFragCounts;
    private final Map<String,Integer> mSampleIndexMap; // map from sampleId into the sample counts matrix

    private static final double ZERO_PREVALENCE_ALLOCATION = 0.03;

    private static final String SAMPLE_ALT_SJ_FILE = "sample_alt_sj_matrix_file";
    private static final String FRAG_COUNT_LOG_VALUE = "alt_sj_log_value";
    private static final String COHORT_SAMPLE_DATA_FILE = "alt_sj_cohort_sample_data";

    public AltSjClassifier(final CuppaConfig config, final SampleDataCache sampleDataCache, final CommandLine cmd)
    {
        mConfig = config;
        mSampleDataCache = sampleDataCache;
        mIsValid = true;

        mRefAltSjPrevalence = Maps.newHashMap();
        mRefAltSjStartPositions = Sets.newHashSet();

        mRefAsjIndexMap = Maps.newHashMap();
        mRefCancerTypes = Lists.newArrayList();
        mRefCancerTypeMatrix = null;

        mSampleIndexMap = Maps.newHashMap();
        mSampleFragCounts = null;

        mFragCountLogValue = Double.parseDouble(cmd.getOptionValue(FRAG_COUNT_LOG_VALUE, "-1"));

        if(config.RefAltSjPrevFile.isEmpty())
            return;

        final List<String> ignoreFields = Lists.newArrayList("GeneId", "Chromosome", "PosStart", "PosEnd");
        loadRefAltSjIndices(config.RefAltSjPrevFile);
        mRefCancerTypeMatrix = loadMatrixDataFile(config.RefAltSjPrevFile, mRefCancerTypes, ignoreFields);

        if(mRefCancerTypeMatrix ==  null)
        {
            mIsValid = false;
            return;
        }

        // convert to log values if required
        convertMatrixFragCounts(mRefCancerTypeMatrix);
        mRefCancerTypeMatrix.cacheTranspose();

        // load a cohort sample frag counts matrix if available
        if(cmd.hasOption(SAMPLE_ALT_SJ_FILE))
        {
            mSampleFragCounts = loadMatrixDataFile(cmd.getOptionValue(SAMPLE_ALT_SJ_FILE), mSampleIndexMap, ignoreFields);

            if(mSampleFragCounts ==  null)
            {
                mIsValid = false;
                return;
            }

            convertMatrixFragCounts(mSampleFragCounts);
            mSampleFragCounts.cacheTranspose();

            if(mSampleFragCounts.Rows != mRefCancerTypeMatrix.Rows)
            {
                CUP_LOGGER.error("alt-SJ matrix mismatch between ref({}) and cohort({})",
                        mRefCancerTypeMatrix.Rows, mSampleFragCounts.Rows);

                mIsValid = false;
            }
        }
        else
        {
            mSampleFragCounts = new Matrix(mRefCancerTypeMatrix.Rows, 1);
        }

        // as a prevalence file
        // mIsValid &= loadRefAltSJs(config.RefAltSjPrevFile);
    }

    public static void addCmdLineArgs(Options options)
    {
        options.addOption(SAMPLE_ALT_SJ_FILE, true, "Cohort sample RNA alt-SJ frag counts matrix file");
        options.addOption(FRAG_COUNT_LOG_VALUE, true, "Use log of frag counts plus this value");
    }

    public CategoryType categoryType() { return ALT_SJ; }
    public boolean isValid() { return mIsValid; }
    public void close() {}

    public void processSample(final SampleData sample, final List<SampleResult> results, final List<SampleSimilarity> similarities)
    {
        if(!mIsValid || (mRefAltSjPrevalence.isEmpty() && mRefCancerTypeMatrix == null))
            return;

        if(mSampleIndexMap.isEmpty())
        {
            if(!loadSampleAltSJsToArray(sample.Id))
                return;
        }

        addCancerCssResults(sample, results);

        /*
        final List<AltSjPrevData> sampleAltSJs = loadSampleAltSJs(sample.Id);

        if(sampleAltSJs == null || sampleAltSJs.isEmpty())
            return;

        calcCancerTypeProbability(sample, sampleAltSJs, results);
        */
    }

    // CSS METHOD
    private void addCancerCssResults(final SampleData sample, final List<SampleResult> results)
    {
        int refCancerCount = mRefCancerTypeMatrix.Cols;

        final Map<String,Double> cancerCssTotals = Maps.newHashMap();

        Integer sampleIndex = mSampleFragCounts.Cols > 1 ? mSampleIndexMap.get(sample.Id) : 0;

        if(sampleIndex == null)
        {
            CUP_LOGGER.warn("sample({}) alt-SJ frag counts from matrix not found", sample.Id);
            return;
        }

        final double[] sampleFragCounts = mSampleFragCounts.getCol(sampleIndex);

        for(int i = 0; i < refCancerCount; ++i)
        {
            final String refCancerType = mRefCancerTypes.get(i);

            if(!isKnownCancerType(refCancerType))
                continue;

            if(!checkIsValidCancerType(sample, refCancerType, cancerCssTotals))
                continue;

            boolean matchesCancerType = sample.cancerType().equals(refCancerType);

            final double[] refAsjFragCounts = sample.isRefSample() && matchesCancerType ?
                    adjustRefCounts(mRefCancerTypeMatrix.getCol(i), sampleFragCounts) : mRefCancerTypeMatrix.getCol(i);

            // skip any alt-SJ not present in the sample, but not the reverse
            double css = calcCosineSim(refAsjFragCounts, sampleFragCounts, false, true);

            if(css < RNA_GENE_EXP_CSS_THRESHOLD)
                continue;

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
                        sample.Id, sample.cancerType(), refCancerType, String.format("%.4f", css),
                        sampleSites, refSites, matchedSites);
            }

            double cssWeight = pow(RNA_ALT_SJ_DIFF_EXPONENT, -100 * (1 - css));

            double weightedCss = css * cssWeight;
            cancerCssTotals.put(refCancerType, weightedCss);
        }

        double totalCss = cancerCssTotals.values().stream().mapToDouble(x -> x).sum();

        for(Map.Entry<String,Double> entry : cancerCssTotals.entrySet())
        {
            cancerCssTotals.put(entry.getKey(), entry.getValue() / totalCss);
        }

        results.add(new SampleResult(
                sample.Id, CLASSIFIER, LIKELIHOOD, ALT_SJ.toString(), String.format("%.4g", totalCss), cancerCssTotals));
    }

    private double[] adjustRefCounts(final double[] refCounts, final double[] sampleCounts)
    {
        // remove the sample's counts after first de-logging both counts if required
        double[] adjustedCounts = new double[refCounts.length];

        for(int b = 0; b < refCounts.length; ++b)
        {
            if(mFragCountLogValue < 0)
            {
                adjustedCounts[b] = max(refCounts[b] - sampleCounts[b], 0);
            }
            else if(refCounts[b] == 0 || sampleCounts[b] == 0)
            {
                adjustedCounts[b] = refCounts[b];
            }
            else
            {
                double sampleCount = exp(sampleCounts[b]) - mFragCountLogValue;
                double refCount = exp(refCounts[b]) - mFragCountLogValue;
                adjustedCounts[b] = convertFragCount(max(refCount - sampleCount, 0));
            }
        }

        return adjustedCounts;
    }

    private static String formAltSjKey(final String chromosome, int posStart, int posEnd)
    {
        return String.format("%s-%d-%d", chromosome, posStart, posEnd);
    }

    private double convertFragCount(double fragCount)
    {
        if(mFragCountLogValue < 0)
            return fragCount;

        return fragCount >= 1 ? log(fragCount + mFragCountLogValue) : 0;
    }

    private void convertMatrixFragCounts(final Matrix matrix)
    {
        if(mFragCountLogValue < 0)
            return;

        final double[][] refData = matrix.getData();
        for(int r = 0; r < matrix.Rows; ++r)
        {
            for(int c = 0; c < matrix.Cols; ++c)
            {
                refData[r][c] = convertFragCount(refData[r][c]);
            }
        }
    }

    private boolean loadRefAltSjIndices(final String filename)
    {
        try
        {
            BufferedReader fileReader = new BufferedReader(new FileReader(filename));

            String header = fileReader.readLine();
            final Map<String,Integer> fieldsIndexMap = createFieldsIndexMap(header, DATA_DELIM);

            int chrIndex = fieldsIndexMap.get("Chromosome");
            int posStartIndex = fieldsIndexMap.get("PosStart");
            int posEndIndex = fieldsIndexMap.get("PosEnd");

            String line = fileReader.readLine();
            int altSjIndex = 0;

            while(line != null)
            {
                final String[] items = line.split(DATA_DELIM, -1);

                final String asjKey = formAltSjKey(
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
        final String filename = mConfig.SampleDataDir + sampleId + ALT_SJ_FILE_ID;

        mSampleFragCounts.initialise(0);
        final double[][] sampleFragCounts = mSampleFragCounts.getData();

        if(!Files.exists(Paths.get(filename)))
            return false;

        try
        {
            final List<String> lines = Files.readAllLines(Paths.get(filename));

            final Map<String,Integer> fieldsIndexMap = createFieldsIndexMap(lines.get(0), DATA_DELIM);
            lines.remove(0);

            int chrIndex = fieldsIndexMap.get("Chromosome");
            int geneIdIndex = fieldsIndexMap.get("GeneId");
            int posStartIndex = fieldsIndexMap.get("SjStart");
            int posEndIndex = fieldsIndexMap.get("SjEnd");
            int fragCountIndex = fieldsIndexMap.get("FragCount");

            boolean hasRefAltSJs = false;

            for(String data : lines)
            {
                final String items[] = data.split(DATA_DELIM, -1);

                String chromosome = items[chrIndex];
                String geneId = items[geneIdIndex];
                int posStart = Integer.parseInt(items[posStartIndex]);
                int posEnd = Integer.parseInt(items[posEndIndex]);

                final String asjKey = formAltSjKey(chromosome, posStart, posEnd);

                Integer matrixIndex = mRefAsjIndexMap.get(asjKey);

                if(matrixIndex == null)
                    continue;

                int fragCount = Integer.parseInt(items[fragCountIndex]);
                hasRefAltSJs = true;

                sampleFragCounts[matrixIndex][0] = convertFragCount(fragCount);
            }

            return hasRefAltSJs;
        }
        catch(IOException e)
        {
            CUP_LOGGER.error("failed to load alt splice junction file({}): {}", filename.toString(), e.toString());
            return false;
        }
    }


    // PREVALENCE METHOD

    private AltSjPrevData findRefAltSjData(final AltSjPrevData altSJ)
    {
        final List<AltSjPrevData> refAltSJs = mRefAltSjPrevalence.get(altSJ.Location.Chromosome);

        if(refAltSJs == null)
            return null;

        return refAltSJs.stream().filter(x -> x.matches(altSJ)).findFirst().orElse(null);
    }

    private void calcCancerTypeProbability(
            final SampleData sample, final List<AltSjPrevData> sampleAltSJs, final List<SampleResult> results)
    {
        // taking the set of drivers as a group, report on the combined probability for each cancer type
        final Map<String,Double> cancerProbTotals = Maps.newHashMap();

        final Set<String> cancerTypes = mSampleDataCache.RefCancerSampleData.keySet();
        int matchedAltSJs = 0;

        for(final AltSjPrevData altSJ : sampleAltSJs)
        {
            // check the alt-SJ is one of the reference ones
            final AltSjPrevData refAltSJ = findRefAltSjData(altSJ);

            if(refAltSJ == null)
                continue;

            ++matchedAltSJs;

            for(String cancerType : cancerTypes)
            {
                if(!checkIsValidCancerType(sample, cancerType, cancerProbTotals))
                    continue;

                Double prevalence = refAltSJ.CancerPrevalences.get(cancerType);

                Double probabilityTotal = cancerProbTotals.get(cancerType);

                if(probabilityTotal == null)
                    probabilityTotal = 1.0;

                boolean adjustMatchingCancerPrev = sample.cancerType().equals(cancerType);

                double driverPrevValue;
                double prevalenceTotal = refAltSJ.PrevalenceTotal;

                if(prevalence != null)
                {
                    driverPrevValue = prevalence;

                    if(adjustMatchingCancerPrev)
                    {
                        int cohortSize = mSampleDataCache.getCancerSampleCount(cancerType);
                        double adjustedIncidence = max(driverPrevValue * cohortSize - 1, 0.0);
                        double adjustedDriverPrevValue = cohortSize > 1 ? adjustedIncidence / (cohortSize - 1) : 0;
                        prevalenceTotal -= driverPrevValue - adjustedDriverPrevValue;
                        driverPrevValue = adjustedDriverPrevValue;
                    }
                }
                else
                {
                    driverPrevValue = refAltSJ.MinPrevalence;
                }

                probabilityTotal *= driverPrevValue / prevalenceTotal;
                cancerProbTotals.put(cancerType, probabilityTotal);
            }

            adjustLowProbabilities(cancerProbTotals);
        }

        if(matchedAltSJs > 0)
        {
            final String altSjStr = String.format("ALT-SJs=%d", matchedAltSJs);

            SampleResult result = new SampleResult(
                    sample.Id, ALT_SJ, LIKELIHOOD, ClassifierType.ALT_SJ.toString(), altSjStr, cancerProbTotals);

            results.add(result);
        }
    }

    private boolean loadRefAltSJs(final String filename)
    {
        try
        {
            final List<String> lines = Files.readAllLines(Paths.get(filename));

            final Map<String, Integer> fieldsIndexMap = createFieldsIndexMap(lines.get(0), ",");
            lines.remove(0);

            // CancerType,GeneId,Chromosome,SjStart,SjEnd,Type,Prev
            int cancerIndex = fieldsIndexMap.get("CancerType");
            int geneIdIndex = fieldsIndexMap.get("GeneId");
            int chromosomeIndex = fieldsIndexMap.get("Chromosome");
            int posStartIndex = fieldsIndexMap.get("SjStart");
            int posEndIndex = fieldsIndexMap.get("SjEnd");
            Integer typeIndex = fieldsIndexMap.get("Type"); // for info sake only
            int prevIndex = fieldsIndexMap.get("Prev");

            double noPrevalence = ZERO_PREVALENCE_ALLOCATION / mSampleDataCache.RefCancerSampleData.size();

            for(String data : lines)
            {
                final String items[] = data.split(",", -1);

                String cancerType = items[cancerIndex];
                String chromosome = items[chromosomeIndex];
                String asjType = typeIndex != null ? items[typeIndex] : "N/A";
                int asjPosStart = Integer.parseInt(items[posStartIndex]);
                int asjPosEnd = Integer.parseInt(items[posEndIndex]);
                double prevalence = Double.parseDouble(items[prevIndex]);

                AltSjPrevData altSJ = new AltSjPrevData(items[geneIdIndex], asjType, chromosome, asjPosStart, asjPosEnd, 0);

                AltSjPrevData matchedAltSJ = null;

                List<AltSjPrevData> altSJs = mRefAltSjPrevalence.get(chromosome);

                if(altSJs == null)
                {
                    altSJs = Lists.newArrayList();
                    mRefAltSjPrevalence.put(chromosome, altSJs);
                }
                else
                {
                    matchedAltSJ = altSJs.stream().filter(x -> x.matches(altSJ)).findFirst().orElse(null);
                }

                if(matchedAltSJ == null)
                {
                    altSJs.add(altSJ);
                    matchedAltSJ = altSJ;
                }

                mRefAltSjStartPositions.add(asjPosStart);

                double adjPrevalence = prevalence + noPrevalence;
                matchedAltSJ.CancerPrevalences.put(cancerType, adjPrevalence);
                matchedAltSJ.PrevalenceTotal += adjPrevalence;
                matchedAltSJ.MinPrevalence = noPrevalence;
            }

            // set default value for missing cancer types
            final Set<String> cancerTypes = mSampleDataCache.RefCancerSampleData.keySet();

            for(List<AltSjPrevData> altSJs : mRefAltSjPrevalence.values())
            {
                for(AltSjPrevData altSJ : altSJs)
                {
                    for(String cancerType : cancerTypes)
                    {
                        if(altSJ.CancerPrevalences.containsKey(cancerType))
                            continue;

                        altSJ.PrevalenceTotal += noPrevalence;
                    }
                }
            }

            CUP_LOGGER.info("loaded {} ref alt splice junctions", mRefAltSjPrevalence.values().stream().mapToInt(x -> x.size()).sum());
        }
        catch(IOException e)
        {
            CUP_LOGGER.error("failed to load ref alt splice junction file({}): {}", filename.toString(), e.toString());
            return false;
        }

        return true;
    }

    private static final String ALT_SJ_FILE_ID = ".isf.alt_splice_junc.csv";
    private static final int ALT_SJ_FRAG_COUNT_THRESHOLD = 3;

    private List<AltSjPrevData> loadSampleAltSJs(final String sampleId)
    {
        final String filename = mConfig.SampleDataDir + sampleId + ALT_SJ_FILE_ID;

        final List<AltSjPrevData> sampleAltSJs = Lists.newArrayList();

        if(!Files.exists(Paths.get(filename)))
            return sampleAltSJs;

        try
        {
            final List<String> lines = Files.readAllLines(Paths.get(filename));

            final Map<String,Integer> fieldsIndexMap = createFieldsIndexMap(lines.get(0), DATA_DELIM);
            lines.remove(0);

            int geneIdIndex = fieldsIndexMap.get("GeneId");
            int chromosomeIndex = fieldsIndexMap.get("Chromosome");
            int posStartIndex = fieldsIndexMap.get("SjStart");
            int posEndIndex = fieldsIndexMap.get("SjEnd");
            int typeIndex = fieldsIndexMap.get("Type");
            int fragCountIndex = fieldsIndexMap.get("FragCount");

            for(String data : lines)
            {
                final String items[] = data.split(DATA_DELIM, -1);

                int asjPosStart = Integer.parseInt(items[posStartIndex]);

                if(!mRefAltSjStartPositions.contains(asjPosStart))
                    continue;

                int fragCount = Integer.parseInt(items[fragCountIndex]);

                if(fragCount < ALT_SJ_FRAG_COUNT_THRESHOLD)
                    continue;

                String chromosome = items[chromosomeIndex];
                int asjPosEnd = Integer.parseInt(items[posEndIndex]);

                sampleAltSJs.add(new AltSjPrevData(items[geneIdIndex], items[typeIndex], chromosome, asjPosStart, asjPosEnd));
            }
        }
        catch(IOException e)
        {
            CUP_LOGGER.error("failed to load alt splice junction file({}): {}", filename.toString(), e.toString());
            return null;
        }

        return sampleAltSJs;
    }
}
