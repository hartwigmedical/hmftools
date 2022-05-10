package com.hartwig.hmftools.cup.somatics;

import static com.hartwig.hmftools.common.sigs.DataUtils.convertList;
import static com.hartwig.hmftools.common.stats.Percentiles.PERCENTILE_COUNT;
import static com.hartwig.hmftools.common.stats.Percentiles.buildPercentiles;
import static com.hartwig.hmftools.common.sigs.SnvSigUtils.SNV_TRINUCLEOTIDE_BUCKET_COUNT;
import static com.hartwig.hmftools.common.sigs.SnvSigUtils.populateBucketMap;
import static com.hartwig.hmftools.common.utils.VectorUtils.sumVector;
import static com.hartwig.hmftools.common.utils.FileWriterUtils.closeBufferedWriter;
import static com.hartwig.hmftools.common.utils.FileWriterUtils.createBufferedWriter;
import static com.hartwig.hmftools.common.variant.VariantType.SNP;
import static com.hartwig.hmftools.cup.CuppaConfig.CUP_LOGGER;
import static com.hartwig.hmftools.cup.CuppaConfig.DATA_DELIM;
import static com.hartwig.hmftools.cup.CuppaRefFiles.COHORT_REF_FILE_SIG_DATA_FILE;
import static com.hartwig.hmftools.cup.CuppaRefFiles.REF_FILE_CANCER_POS_FREQ_COUNTS;
import static com.hartwig.hmftools.cup.CuppaRefFiles.REF_FILE_CN_ADJUSTED_SAMPLE_POS_FREQ_COUNTS;
import static com.hartwig.hmftools.cup.CuppaRefFiles.REF_FILE_COPY_NUMBER_PROFILE;
import static com.hartwig.hmftools.cup.CuppaRefFiles.REF_FILE_SAMPLE_POS_FREQ_COUNTS;
import static com.hartwig.hmftools.cup.CuppaRefFiles.REF_FILE_SIG_PERC;
import static com.hartwig.hmftools.cup.CuppaRefFiles.REF_FILE_SNV_COUNTS;
import static com.hartwig.hmftools.cup.CuppaRefFiles.purpleSomaticVcfFile;
import static com.hartwig.hmftools.cup.common.CategoryType.SNV;
import static com.hartwig.hmftools.cup.common.CupConstants.CANCER_TYPE_OTHER;
import static com.hartwig.hmftools.cup.common.CupConstants.POS_FREQ_BUCKET_SIZE;
import static com.hartwig.hmftools.cup.common.CupConstants.POS_FREQ_MAX_SAMPLE_COUNT;
import static com.hartwig.hmftools.cup.common.SampleData.isKnownCancerType;
import static com.hartwig.hmftools.cup.ref.RefDataConfig.parseFileSet;
import static com.hartwig.hmftools.cup.somatics.CopyNumberProfile.buildCopyNumberProfile;
import static com.hartwig.hmftools.cup.somatics.GenomicPositions.addPanCancerNoise;
import static com.hartwig.hmftools.cup.somatics.GenomicPositions.buildCancerMatrix;
import static com.hartwig.hmftools.cup.somatics.GenomicPositions.extractPositionFrequencyCounts;
import static com.hartwig.hmftools.cup.somatics.SomaticClassifier.SNV_POS_FREQ_POS_SIZE;
import static com.hartwig.hmftools.cup.somatics.SomaticDataLoader.loadRefSampleCounts;
import static com.hartwig.hmftools.cup.somatics.SomaticDataLoader.loadSigContribsFromCohortFile;
import static com.hartwig.hmftools.cup.somatics.SomaticDataLoader.loadSomaticVariants;
import static com.hartwig.hmftools.cup.somatics.SomaticSigs.SIG_NAME_13;
import static com.hartwig.hmftools.cup.somatics.SomaticSigs.SIG_NAME_2;
import static com.hartwig.hmftools.cup.somatics.SomaticsCommon.DEC_3_FORMAT;
import static com.hartwig.hmftools.cup.somatics.SomaticsCommon.INTEGER_FORMAT;
import static com.hartwig.hmftools.cup.somatics.SomaticsCommon.NORMALISE_COPY_NUMBER;
import static com.hartwig.hmftools.cup.somatics.SomaticsCommon.NORMALISE_COPY_NUMBER_DESC;
import static com.hartwig.hmftools.cup.somatics.SomaticsCommon.EXCLUDE_AID_APOBEC;
import static com.hartwig.hmftools.cup.somatics.SomaticsCommon.EXCLUDE_AID_APOBEC_DESC;
import static com.hartwig.hmftools.cup.somatics.SomaticsCommon.loadMultipleMatrixFiles;
import static com.hartwig.hmftools.cup.somatics.SomaticsCommon.writeMatrix;
import static com.hartwig.hmftools.cup.somatics.SomaticsCommon.writeSampleMatrix;
import static com.hartwig.hmftools.cup.somatics.TrinucleotideCounts.extractTrinucleotideCounts;

import java.io.BufferedWriter;
import java.io.File;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.List;
import java.util.Map;
import java.util.stream.Collectors;

import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import com.hartwig.hmftools.common.sigs.PositionFrequencies;
import com.hartwig.hmftools.common.utils.Matrix;
import com.hartwig.hmftools.cup.common.CategoryType;
import com.hartwig.hmftools.cup.common.SampleData;
import com.hartwig.hmftools.cup.common.SampleDataCache;
import com.hartwig.hmftools.cup.ref.RefDataConfig;
import com.hartwig.hmftools.cup.ref.RefClassifier;
import com.hartwig.hmftools.cup.traits.SampleTraitsData;

import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.Options;
import org.jetbrains.annotations.NotNull;

public class RefSomatics implements RefClassifier
{
    private final RefDataConfig mConfig;
    private final SampleDataCache mSampleDataCache;

    private final Map<String,Map<String,List<Double>>> mCancerSigContribs;
    private final Map<String,List<Double>> mCancerSnvCounts;

    private Matrix mTriNucCounts; // counts per tri-nucleotide bucket
    private final Map<String,Integer> mTriNucCountsIndex;
    private boolean mWriteTriNucMatrixData;

    private Matrix mPosFreqCounts; // counts per genomic position
    private final Map<String,Integer> mPosFreqCountsIndex;
    private boolean mWritePosFreqMatrixData;

    private final boolean mExcludeAidApobec;

    private final boolean mBuildCopyNumber;
    private final boolean mApplyCopyNumber; // to genomic positions
    private Matrix mCopyNumberProfile; // same dimensions as genomic position
    private final int mGenPosNoiseAllocation;

    private final PositionFrequencies mPosFrequencies;

    private BufferedWriter mRefDataWriter;

    public static final String REF_SIG_TYPE_SNV_COUNT = "SnvCount";
    private static final String MATRIX_TYPE_SNV_96 = "SNV_96";
    private static final String MATRIX_TYPE_GEN_POS = "GEN_POS";
    private static final String MATRIX_TYPE_CN_PROFILE = "CN_PROFILE";

    // config
    private static final String BUILD_CN_PROFILE = "build_copy_number";
    public static final String GEN_POS_NOISE_ALLOC = "gen_pos_noise_alloc";
    public static final String GEN_POS_NOISE_ALLOC_DESC = "Add genomic position counts from pan-cancer medians";

    public RefSomatics(final RefDataConfig config, final SampleDataCache sampleDataCache, final CommandLine cmd)
    {
        mConfig = config;
        mSampleDataCache = sampleDataCache;

        mCancerSigContribs = Maps.newHashMap();
        mCancerSnvCounts = Maps.newHashMap();

        mRefDataWriter = null;

        mTriNucCounts = null;
        mTriNucCountsIndex = Maps.newHashMap();
        mWriteTriNucMatrixData = false;

        mPosFreqCounts = null;
        mPosFreqCountsIndex = Maps.newHashMap();
        mWritePosFreqMatrixData = false;

        mExcludeAidApobec = cmd.hasOption(EXCLUDE_AID_APOBEC);
        mApplyCopyNumber = cmd.hasOption(NORMALISE_COPY_NUMBER);

        int posFreqBucketSize = cmd.hasOption(SNV_POS_FREQ_POS_SIZE) ?
                Integer.parseInt(cmd.getOptionValue(SNV_POS_FREQ_POS_SIZE)) : POS_FREQ_BUCKET_SIZE;

        mPosFrequencies = new PositionFrequencies(posFreqBucketSize, POS_FREQ_MAX_SAMPLE_COUNT);

        mBuildCopyNumber = cmd.hasOption(BUILD_CN_PROFILE);
        mGenPosNoiseAllocation = Integer.parseInt(cmd.getOptionValue(GEN_POS_NOISE_ALLOC, "0"));
        mCopyNumberProfile = null;
    }

    public CategoryType categoryType() { return SNV; }

    public static boolean requiresBuild(final RefDataConfig config)
    {
        if(config.DbAccess != null)
            return true;

        if(config.Categories.contains(SNV))
            return true;

        return !config.CohortSigContribsFile.isEmpty() || !config.Snv96MatrixFile.isEmpty() || !config.PurpleDir.isEmpty();
    }

    public static void addCmdLineArgs(@NotNull Options options)
    {
        options.addOption(SNV_POS_FREQ_POS_SIZE, true, "Genomic position bucket size (default: 20000)");
        options.addOption(EXCLUDE_AID_APOBEC, false, EXCLUDE_AID_APOBEC_DESC);
        options.addOption(NORMALISE_COPY_NUMBER, false, NORMALISE_COPY_NUMBER_DESC);
        options.addOption(BUILD_CN_PROFILE, false, "Build a copy-number profile to match genomic positions");
        options.addOption(GEN_POS_NOISE_ALLOC, true, GEN_POS_NOISE_ALLOC_DESC);
    }

    public void buildRefDataSets()
    {
        CUP_LOGGER.info("building SNV and signatures reference data");

        final List<String> snvCountFiles = parseFileSet(mConfig.Snv96MatrixFile);
        final List<String> posFreqCountFiles = parseFileSet(mConfig.GenPosMatrixFile);
        final List<String> cnProfileFiles = parseFileSet(mConfig.CopyNumberMatrixFile);

        boolean combineCohortFiles = snvCountFiles.size() > 1 && posFreqCountFiles.size() == snvCountFiles.size();

        if(combineCohortFiles)
        {
            final List<String> refSampleIds = mSampleDataCache.refSampleIds(false);
            mTriNucCounts = loadMultipleMatrixFiles(snvCountFiles, refSampleIds, mTriNucCountsIndex, MATRIX_TYPE_SNV_96);
            mPosFreqCounts = loadMultipleMatrixFiles(posFreqCountFiles, refSampleIds, mPosFreqCountsIndex, MATRIX_TYPE_GEN_POS);
            mCopyNumberProfile = loadMultipleMatrixFiles(cnProfileFiles, refSampleIds, mPosFreqCountsIndex, MATRIX_TYPE_CN_PROFILE);
            mWriteTriNucMatrixData = mWritePosFreqMatrixData = true;
        }
        else
        {
            mTriNucCounts = loadReferenceSnvCounts(mConfig.Snv96MatrixFile, mTriNucCountsIndex, MATRIX_TYPE_SNV_96);
            mPosFreqCounts = loadReferenceSnvCounts(mConfig.GenPosMatrixFile, mPosFreqCountsIndex, MATRIX_TYPE_GEN_POS);

            retrieveMissingSampleCounts();

            if(!mConfig.CopyNumberMatrixFile.isEmpty())
            {
                final List<String> existingRefSampleIds = Lists.newArrayList();
                mCopyNumberProfile = loadRefSampleCounts(mConfig.CopyNumberMatrixFile, existingRefSampleIds, Lists.newArrayList());
            }
            else if(mBuildCopyNumber)
            {
                final List<String> refSampleIds = mSampleDataCache.refSampleIds(false);
                mCopyNumberProfile = buildCopyNumberProfile(refSampleIds, mConfig, mPosFreqCountsIndex, mPosFrequencies);
            }
        }

        mTriNucCounts.cacheTranspose();
        mPosFreqCounts.cacheTranspose();

        // write out sample matrix data unless they were already correct
        if(mWriteTriNucMatrixData)
            writeSampleMatrix(mTriNucCounts, mTriNucCountsIndex, mConfig.OutputDir + REF_FILE_SNV_COUNTS, INTEGER_FORMAT);

        if(mWritePosFreqMatrixData)
        {
            writeSampleMatrix(mPosFreqCounts, mPosFreqCountsIndex, mConfig.OutputDir + REF_FILE_SAMPLE_POS_FREQ_COUNTS, INTEGER_FORMAT);
        }

        if(mBuildCopyNumber)
        {
            writeSampleMatrix(mCopyNumberProfile, mPosFreqCountsIndex, mConfig.OutputDir + REF_FILE_COPY_NUMBER_PROFILE, DEC_3_FORMAT);
        }

        Matrix sampleGenPosCounts = null;

        if(mApplyCopyNumber)
        {
            CUP_LOGGER.info("normalising sample pos-freq counts with copy-number");

            Matrix sampleCnNormalisedCounts = CopyNumberProfile.buildSampleCopyNumberNormalisedCounts(
                    mPosFreqCounts, mPosFreqCountsIndex, mCopyNumberProfile, mSampleDataCache.RefSampleDataList,
                    mSampleDataCache.RefSampleTraitsData);

            sampleGenPosCounts = sampleCnNormalisedCounts;

            if(mBuildCopyNumber)
            {
                writeSampleMatrix(
                        sampleCnNormalisedCounts, mPosFreqCountsIndex,
                        mConfig.OutputDir + REF_FILE_CN_ADJUSTED_SAMPLE_POS_FREQ_COUNTS, DEC_3_FORMAT);
            }
        }
        else
        {
            sampleGenPosCounts = mPosFreqCounts;
        }

        buildSignaturePercentiles();
        buildSnvCountPercentiles();

        List<String> cancerTypes = mSampleDataCache.RefCancerSampleData.keySet().stream()
                .filter(x -> !x.equals(CANCER_TYPE_OTHER)).collect(Collectors.toList());

        Matrix cancerGenPosMatrix = buildCancerMatrix(
                sampleGenPosCounts, mPosFreqCountsIndex, cancerTypes, mSampleDataCache.RefCancerSampleData, mPosFrequencies.getMaxSampleCount());

        if(mGenPosNoiseAllocation > 0)
        {
            addPanCancerNoise(cancerGenPosMatrix, mGenPosNoiseAllocation);
        }

        writeMatrix(cancerGenPosMatrix, cancerTypes, mConfig.OutputDir + REF_FILE_CANCER_POS_FREQ_COUNTS, DEC_3_FORMAT);

        closeBufferedWriter(mRefDataWriter);
    }

    private Matrix loadReferenceSnvCounts(final String refFilename, final Map<String,Integer> sampleCountsIndex, final String type)
    {
        if(refFilename.isEmpty())
        {
            if(type.equals(MATRIX_TYPE_SNV_96))
                mWriteTriNucMatrixData = true;
            else
                mWritePosFreqMatrixData = true;

            return null;
        }

        CUP_LOGGER.debug("loading SNV {} reference data", type);

        // check if complete file has already been provided (eg if only other reference data is being built)
        Matrix refMatrix = null;

        final List<String> existingRefSampleIds = Lists.newArrayList();
        final Matrix existingRefSampleCounts = loadRefSampleCounts(refFilename, existingRefSampleIds, Lists.newArrayList("BucketName"));

        final List<String> refSampleIds = mSampleDataCache.refSampleIds(false);
        boolean hasMissingSamples = refSampleIds.stream().anyMatch(x -> !existingRefSampleIds.contains(x));

        if(!hasMissingSamples && existingRefSampleIds.size() == refSampleIds.size())
        {
            CUP_LOGGER.debug("using existing SNV {} reference counts", type);

            for(int i = 0; i < existingRefSampleIds.size(); ++i)
            {
                sampleCountsIndex.put(existingRefSampleIds.get(i), i);
            }

            if(type.equals(MATRIX_TYPE_SNV_96))
                mWriteTriNucMatrixData = false;
            else
                mWritePosFreqMatrixData = false;

            return existingRefSampleCounts;
        }

        // take any existing counts
        if(existingRefSampleCounts != null)
        {
            existingRefSampleCounts.cacheTranspose();

            refMatrix = new Matrix(existingRefSampleCounts.Rows, refSampleIds.size());

            int refSampleIndex = 0;

            for(int i = 0; i < existingRefSampleIds.size(); ++i)
            {
                final String sampleId = existingRefSampleIds.get(i);

                if(!mSampleDataCache.hasRefSample(sampleId))
                    continue;

                refMatrix.setCol(refSampleIndex, existingRefSampleCounts.getCol(i));
                sampleCountsIndex.put(existingRefSampleIds.get(i), refSampleIndex);
                ++refSampleIndex;
            }
        }

        if(type.equals(MATRIX_TYPE_SNV_96))
            mWriteTriNucMatrixData = true;
        else
            mWritePosFreqMatrixData = true;

        return refMatrix;
    }

    private void retrieveMissingSampleCounts()
    {
        // returns true if both types of counts are already exactly accounted for in the loaded matrix data
        final List<String> refSampleIds = mSampleDataCache.refSampleIds(false);

        long missingSamples = refSampleIds.stream()
                .filter(x -> !mTriNucCountsIndex.containsKey(x) || !mPosFreqCountsIndex.containsKey(x)).count();

        if(missingSamples == 0)
            return;

        int refSampleCount = refSampleIds.size();
        CUP_LOGGER.debug("retrieving SNV data for {} samples from refSampleCount({})", missingSamples, refSampleCount);

        if(mTriNucCounts == null)
        {
            mTriNucCounts = new Matrix(SNV_TRINUCLEOTIDE_BUCKET_COUNT, refSampleCount);
        }

        if(mPosFreqCounts == null)
        {
            mPosFreqCounts = new Matrix(mPosFrequencies.getBucketCount(), refSampleCount);
        }

        final Map<String,Integer> triNucBucketNameMap = Maps.newHashMap();
        populateBucketMap(triNucBucketNameMap);

        int nextLog = 100;
        int retrievedSamples = 0;

        for(int i = 0; i < refSampleCount; ++i)
        {
            final String sampleId = refSampleIds.get(i);
            boolean needsTriNucCounts = !mTriNucCountsIndex.containsKey(sampleId);
            boolean needsPosFreqCounts = !mPosFreqCountsIndex.containsKey(sampleId);

            if(!needsPosFreqCounts && !needsTriNucCounts)
                continue;

            ++retrievedSamples;

            if(retrievedSamples >= nextLog)
            {
                nextLog += 100;
                CUP_LOGGER.debug("retrieved SNV data for {} samples", retrievedSamples);
            }

            List<SomaticVariant> variants = Lists.newArrayList();

            if(mConfig.DbAccess != null)
            {
                variants.addAll(loadSomaticVariants(sampleId, mConfig.DbAccess));
            }
            else
            {
                final String somaticVcfFile = purpleSomaticVcfFile(mConfig.PurpleDir, sampleId);
                variants.addAll(loadSomaticVariants(somaticVcfFile, Lists.newArrayList(SNP)));
            }

            if(needsTriNucCounts)
            {
                final double[] triNucCounts = extractTrinucleotideCounts(variants, triNucBucketNameMap);

                int refSampleIndex = mTriNucCountsIndex.size();
                mTriNucCounts.setCol(refSampleIndex, triNucCounts);
                mTriNucCountsIndex.put(sampleId, refSampleIndex);
            }

            if(needsPosFreqCounts)
            {
                int refSampleIndex = mPosFreqCountsIndex.size();
                mPosFreqCountsIndex.put(sampleId, refSampleIndex);

                AidApobecStatus aidApobecStatus = mExcludeAidApobec ? AidApobecStatus.FALSE_ONLY : AidApobecStatus.ALL;
                extractPositionFrequencyCounts(variants, mPosFrequencies, aidApobecStatus);
                mPosFreqCounts.setCol(refSampleIndex, mPosFrequencies.getCounts());
            }
        }
    }

    private void buildSignaturePercentiles()
    {
        CUP_LOGGER.debug("building signature allocation reference data");

        initialiseRefDataWriter();

        final Map<String,Map<String,Double>> sampleSigContributions = Maps.newHashMap();

        if(!mConfig.CohortSigContribsFile.isEmpty())
        {
            final Map<String,Map<String,Double>> allSampleSigContributions = Maps.newHashMap();

            final List<String> files = parseFileSet(mConfig.CohortSigContribsFile);
            files.forEach(x -> loadSigContribsFromCohortFile(x, allSampleSigContributions));

            // extract only reference sample data
            for(Map.Entry<String,Map<String,Double>> entry : allSampleSigContributions.entrySet())
            {
                String sampleId = entry.getKey();

                if(mSampleDataCache.hasRefSample(sampleId))
                {
                    sampleSigContributions.put(sampleId, entry.getValue());
                }
            }
        }
        else if(mConfig.DbAccess != null)
        {
            SomaticDataLoader.loadSigContribsFromDatabase(
                    mConfig.DbAccess, mSampleDataCache.refSampleIds(true), sampleSigContributions);
        }

        writeCohortData(sampleSigContributions);

        for(Map.Entry<String,Map<String,Double>> entry : sampleSigContributions.entrySet())
        {
            final String sampleId = entry.getKey();
            final Map<String,Double> sigAllocations = entry.getValue();

            final String cancerType = mSampleDataCache.RefSampleCancerTypeMap.get(sampleId);

            if(cancerType == null)
            {
                // expected if a smaller ref sample set is being run
                // CUP_LOGGER.debug("sample({}) signatures missing cancer type", sampleId);
                continue;
            }

            Map<String,List<Double>> sigDataMap = mCancerSigContribs.get(cancerType);

            if(sigDataMap == null)
            {
                sigDataMap = Maps.newHashMap();
                mCancerSigContribs.put(cancerType, sigDataMap);
            }

            for(String sigName : SomaticSigs.REPORTABLE_SIGS.keySet())
            {
                double sigContrib = sigAllocations.containsKey(sigName) ? sigAllocations.get(sigName) : 0;

                // combine 2 & 13
                if(sigName.equalsIgnoreCase(SIG_NAME_13))
                    continue;

                if(sigName.equalsIgnoreCase(SIG_NAME_2))
                {
                    sigContrib += sigAllocations.containsKey(SIG_NAME_13) ? sigAllocations.get(SIG_NAME_13) : 0;
                }

                List<Double> sigContribs = sigDataMap.get(sigName);

                if(sigContribs == null)
                {
                    sigDataMap.put(sigName, Lists.newArrayList(sigContrib));
                    continue;
                }

                // add in ascending order
                int index = 0;
                while(index < sigContribs.size())
                {
                    if(sigContrib < sigContribs.get(index))
                        break;

                    ++index;
                }

                sigContribs.add(index, sigContrib);
            }
        }

        for(Map.Entry<String,Map<String,List<Double>>> entry : mCancerSigContribs.entrySet())
        {
            final String cancerType = entry.getKey();

            if(!mSampleDataCache.hasRefCancerType(cancerType))
                continue;

            for(Map.Entry<String,List<Double>> sigEntry : entry.getValue().entrySet())
            {
                final String sigName = sigEntry.getKey();
                final double[] percentiles = buildPercentiles(convertList(sigEntry.getValue()));
                writeRefSigData(cancerType, sigName, percentiles);
            }
        }
    }

    public static String convertSignatureName(final String sigName)
    {
        return sigName.replaceAll("Signature.", "Sig");
    }

    private void buildSnvCountPercentiles()
    {
        for(Map.Entry<String,Integer> entry : mTriNucCountsIndex.entrySet())
        {
            final String sampleId = entry.getKey();
            double sampleTotal = sumVector(mTriNucCounts.getCol(entry.getValue()));

            final String cancerType = mSampleDataCache.RefSampleCancerTypeMap.get(sampleId);

            if(cancerType == null)
            {
                // not a ref sample even though in the counts file
                CUP_LOGGER.debug("sample({}) SNV missing cancer type", sampleId);
                continue;
            }

            List<Double> sampleCounts = mCancerSnvCounts.get(cancerType);
            if(sampleCounts == null)
            {
                mCancerSnvCounts.put(cancerType, Lists.newArrayList(sampleTotal));
            }
            else
            {
                int index = 0;
                while(index < sampleCounts.size())
                {
                    if(sampleTotal < sampleCounts.get(index))
                        break;

                    ++index;
                }

                sampleCounts.add(index, sampleTotal);
            }
        }

        for(Map.Entry<String,List<Double>> entry : mCancerSnvCounts.entrySet())
        {
            final String cancerType = entry.getKey();

            if(!isKnownCancerType(cancerType))
                continue;

            final double[] percentiles = buildPercentiles(convertList(entry.getValue()));
            writeRefSnvCountData(cancerType, percentiles);
        }
    }
    
    private void initialiseRefDataWriter()
    {
        try
        {
            final String filename = mConfig.OutputDir + REF_FILE_SIG_PERC;
            mRefDataWriter = createBufferedWriter(filename, false);

            mRefDataWriter.write("CancerType,DataType");

            for(int i = 0; i < PERCENTILE_COUNT; ++i)
            {
                mRefDataWriter.write(String.format(",Pct_%.2f", i * 0.01));
            }

            mRefDataWriter.newLine();
        }
        catch(IOException e)
        {
            CUP_LOGGER.error("failed to write signatures ref data output: {}", e.toString());
        }
    }

    private void writeRefSigData(final String cancerType, final String sigName, final double[] percentileValues)
    {
        try
        {
            mRefDataWriter.write(String.format("%s,%s", cancerType, sigName));

            for(int i = 0; i < percentileValues.length; ++i)
            {
                mRefDataWriter.write(String.format(",%.6f", percentileValues[i]));
            }

            mRefDataWriter.newLine();
        }
        catch(IOException e)
        {
            CUP_LOGGER.error("failed to write signatures ref data output: {}", e.toString());
        }
    }

    private void writeRefSnvCountData(final String cancerType, final double[] percentileValues)
    {
        try
        {
            mRefDataWriter.write(String.format("%s,%s", cancerType, REF_SIG_TYPE_SNV_COUNT));

            for(int i = 0; i < percentileValues.length; ++i)
            {
                mRefDataWriter.write(String.format(",%.0f", percentileValues[i]));
            }

            mRefDataWriter.newLine();
        }
        catch(IOException e)
        {
            CUP_LOGGER.error("failed to write signatures ref data output: {}", e.toString());
        }
    }

    public static boolean populateRefPercentileData(
            final String filename, final Map<String,Map<String,double[]>> cancerSigContribs, final Map<String,double[]> cancerSnvCounts)
    {
        try
        {
            final List<String> fileData = Files.readAllLines(new File(filename).toPath());

            final String header = fileData.get(0);
            fileData.remove(0);

            for(final String line : fileData)
            {
                // SampleId,DataType,Pct_0.00 etc
                final String[] items = line.split(DATA_DELIM, -1);
                String cancerType = items[0];

                String dataType = items[1];

                double[] percentileData = new double[PERCENTILE_COUNT];

                int startIndex = 2;

                for(int i = startIndex; i < items.length; ++i)
                {
                    double value = Double.parseDouble(items[i]);
                    percentileData[i - startIndex] = value;
                }

                if(dataType.equals(REF_SIG_TYPE_SNV_COUNT))
                {
                    cancerSnvCounts.put(cancerType, percentileData);
                }
                else
                {
                    String sigName = dataType;

                    Map<String, double[]> sigContribsMap = cancerSigContribs.get(cancerType);

                    if(sigContribsMap == null)
                    {
                        sigContribsMap = Maps.newHashMap();
                        cancerSigContribs.put(cancerType, sigContribsMap);
                    }

                    sigContribsMap.put(sigName, percentileData);
                }
            }
        }
        catch (IOException e)
        {
            CUP_LOGGER.error("failed to read sig contrib percentile data file({}): {}", filename, e.toString());
            return false;
        }

        return true;
    }

    private void writeCohortData(final Map<String,Map<String,Double>> sampleSigContributions)
    {
        if(!mConfig.WriteCohortFiles)
            return;

        final String filename = mConfig.OutputDir + COHORT_REF_FILE_SIG_DATA_FILE;
        if(Files.exists(Paths.get(filename)))
        {
            CUP_LOGGER.warn("not over-writing cohort sig-contributions reference file({})", filename);
            return;
        }

        CUP_LOGGER.info("writing cohort signature allocation reference data");

        try
        {
            BufferedWriter writer = createBufferedWriter(filename, false);

            writer.write(SampleTraitsData.header());
            writer.newLine();

            for(Map.Entry<String,Map<String,Double>> entry : sampleSigContributions.entrySet())
            {
                final String sampleId = entry.getKey();
                final Map<String,Double> sigAllocs = entry.getValue();

                for(Map.Entry<String,Double> sigEntry : sigAllocs.entrySet())
                {
                    final String sigName = sigEntry.getKey();

                    if(SomaticSigs.REPORTABLE_SIGS.keySet().contains(sigName))
                    {
                        writer.write(String.format("%s,%s,%s", sampleId, sigName, sigEntry.getValue()));
                        writer.newLine();
                    }
                }
            }

            closeBufferedWriter(writer);
        }
        catch(IOException e)
        {
            CUP_LOGGER.error("failed to write signature allocation cohort data output: {}", e.toString());
        }
    }
}
