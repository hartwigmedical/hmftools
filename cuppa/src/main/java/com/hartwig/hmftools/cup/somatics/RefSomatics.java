package com.hartwig.hmftools.cup.somatics;

import static com.hartwig.hmftools.common.sigs.DataUtils.convertList;
import static com.hartwig.hmftools.common.sigs.Percentiles.PERCENTILE_COUNT;
import static com.hartwig.hmftools.common.sigs.Percentiles.buildPercentiles;
import static com.hartwig.hmftools.common.sigs.VectorUtils.sumVector;
import static com.hartwig.hmftools.common.utils.io.FileWriterUtils.closeBufferedWriter;
import static com.hartwig.hmftools.common.utils.io.FileWriterUtils.createBufferedWriter;
import static com.hartwig.hmftools.cup.CuppaConfig.CUP_LOGGER;
import static com.hartwig.hmftools.cup.CuppaConfig.DATA_DELIM;
import static com.hartwig.hmftools.cup.common.CategoryType.GENE_EXP;
import static com.hartwig.hmftools.cup.common.CategoryType.SNV;
import static com.hartwig.hmftools.cup.somatics.SomaticDataLoader.loadRefSampleCounts;

import java.io.BufferedWriter;
import java.io.File;
import java.io.IOException;
import java.nio.file.Files;
import java.util.List;
import java.util.Map;

import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import com.hartwig.hmftools.common.sigs.SigMatrix;
import com.hartwig.hmftools.cup.common.CategoryType;
import com.hartwig.hmftools.cup.common.SampleDataCache;
import com.hartwig.hmftools.cup.ref.RefDataConfig;
import com.hartwig.hmftools.cup.rna.RefClassifier;

public class RefSomatics implements RefClassifier
{
    private final RefDataConfig mConfig;
    private final SampleDataCache mSampleDataCache;

    private final Map<String,Map<String,List<Double>>> mCancerSigContribs;
    private final Map<String,List<Double>> mCancerSnvCounts;

    private final List<String> mSampleNames;
    private final SigMatrix mSampleCounts;

    private BufferedWriter mRefDataWriter;

    public static final String REF_SIG_TYPE_CONTRIB = "SigContrib";
    public static final String REF_SIG_TYPE_SNV_COUNT = "SnvCount";

    public RefSomatics(final RefDataConfig config, final SampleDataCache sampleDataCache)
    {
        mConfig = config;
        mSampleDataCache = sampleDataCache;

        mCancerSigContribs = Maps.newHashMap();
        mCancerSnvCounts = Maps.newHashMap();

        mRefDataWriter = null;

        loadRefSigContributions(mConfig.RefSigContribsFile);
        mSampleNames = Lists.newArrayList();
        mSampleCounts = loadRefSampleCounts(mConfig.RefSnvCountsFile, mSampleNames);
    }

    public CategoryType categoryType() { return SNV; }

    public static boolean requiresBuild(final RefDataConfig config)
    {
        return !config.RefSigContribsFile.isEmpty() && !config.RefSnvCountsFile.isEmpty();
    }

    public void buildRefDataSets()
    {
        if(mConfig.RefSigContribsFile.isEmpty() && mConfig.RefSnvCountsFile.isEmpty())
            return;

        CUP_LOGGER.info("building SNV and signatures reference data");

        initialiseRefDataWriter();

        for(Map.Entry<String,Map<String,List<Double>>> entry : mCancerSigContribs.entrySet())
        {
            final String cancerType = entry.getKey();

            for(Map.Entry<String,List<Double>> sigEntry : entry.getValue().entrySet())
            {
                final String sigName = sigEntry.getKey();
                final double[] percentiles = buildPercentiles(convertList(sigEntry.getValue()));
                writeRefSigData(cancerType, sigName, percentiles);
            }
        }

        for(int i = 0; i < mSampleNames.size(); ++i)
        {
            double sampleTotal = sumVector(mSampleCounts.getCol(i));
            final String sampleId = mSampleNames.get(i);
            final String cancerType = mSampleDataCache.RefSampleCancerTypeMap.get(sampleId);

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

            final double[] percentiles = buildPercentiles(convertList(entry.getValue()));
            writeRefSnvCountData(cancerType, percentiles);
        }

        closeBufferedWriter(mRefDataWriter);
    }
    
    private void initialiseRefDataWriter()
    {
        try
        {
            final String filename = mConfig.OutputDir + "cup_ref_sig_percentiles.csv";
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

    private void loadRefSigContributions(final String filename)
    {
        if(filename.isEmpty())
            return;

        try
        {
            final List<String> fileData = Files.readAllLines(new File(filename).toPath());

            final String header = fileData.get(0);
            fileData.remove(0);

            for(final String line : fileData)
            {
                // SampleId,SigName,SigContrib,SigPercent
                final String[] items = line.split(DATA_DELIM, -1);
                String sampleId = items[0];
                String sigName = items[1];
                double sigContrib = Double.parseDouble(items[2]);

                final String cancerType = mSampleDataCache.RefSampleCancerTypeMap.get(sampleId);

                if(cancerType == null)
                {
                    CUP_LOGGER.error("sample({}) signatures missing cancer type", sampleId);
                    continue;
                }

                Map<String,List<Double>> sigDataMap = mCancerSigContribs.get(cancerType);

                if(sigDataMap == null)
                {
                    sigDataMap = Maps.newHashMap();
                    mCancerSigContribs.put(cancerType, sigDataMap);
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
        catch (IOException e)
        {
            CUP_LOGGER.error("failed to read sig contribution data file({}): {}", filename, e.toString());
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
}
