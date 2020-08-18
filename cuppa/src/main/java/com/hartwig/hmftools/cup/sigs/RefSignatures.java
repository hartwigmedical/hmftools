package com.hartwig.hmftools.cup.sigs;

import static com.hartwig.hmftools.common.sigs.DataUtils.convertList;
import static com.hartwig.hmftools.common.sigs.Percentiles.PERCENTILE_COUNT;
import static com.hartwig.hmftools.common.sigs.Percentiles.buildPercentiles;
import static com.hartwig.hmftools.common.utils.io.FileWriterUtils.closeBufferedWriter;
import static com.hartwig.hmftools.common.utils.io.FileWriterUtils.createBufferedWriter;
import static com.hartwig.hmftools.cup.SampleAnalyserConfig.CUP_LOGGER;
import static com.hartwig.hmftools.cup.SampleAnalyserConfig.DATA_DELIM;

import java.io.BufferedWriter;
import java.io.File;
import java.io.IOException;
import java.nio.file.Files;
import java.util.List;
import java.util.Map;

import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import com.hartwig.hmftools.cup.common.SampleDataCache;
import com.hartwig.hmftools.cup.ref.RefDataConfig;

public class RefSignatures
{
    private final RefDataConfig mConfig;
    private final SampleDataCache mSampleDataCache;

    private final Map<String,Map<String,List<Double>>> mCancerSigContribs;

    private BufferedWriter mRefDataWriter;

    public RefSignatures(final RefDataConfig config, final SampleDataCache sampleDataCache)
    {
        mConfig = config;
        mSampleDataCache = sampleDataCache;

        mCancerSigContribs = Maps.newHashMap();
        mRefDataWriter = null;

        initialiseRefDataWriter();

        loadRefSigContributions(mConfig.RefSigContribsFile);
    }

    public void buildRefDataSets()
    {
        buildRefSignatureData();
    }

    private void buildRefSignatureData()
    {
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

        closeBufferedWriter(mRefDataWriter);
    }
    
    private void initialiseRefDataWriter()
    {
        try
        {
            final String filename = mConfig.OutputDir + "cup_ref_sig_percentiles.csv";
            mRefDataWriter = createBufferedWriter(filename, false);

            mRefDataWriter.write("CancerType,SigName");

            for(int i = 0; i < PERCENTILE_COUNT; ++i)
            {
                mRefDataWriter.write(String.format(",Pct_%.2f", i * 0.01));
            }

            mRefDataWriter.newLine();
        }
        catch(IOException e)
        {
            CUP_LOGGER.error("failed to write sample traits ref data output: {}", e.toString());
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
            CUP_LOGGER.error("failed to write sample traits ref data output: {}", e.toString());
        }
    }

    private void loadRefSigContributions(final String filename)
    {
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
                    CUP_LOGGER.error("sample({}) missing cancer type", sampleId);
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
            CUP_LOGGER.error("failed to read driver prevalence data file({}): {}", filename, e.toString());
        }
    }

    public static void populateRefSigContributions(final String filename, final Map<String,Map<String,double[]>> cancerSigContribsMap)
    {
        try
        {
            final List<String> fileData = Files.readAllLines(new File(filename).toPath());

            final String header = fileData.get(0);
            fileData.remove(0);

            for(final String line : fileData)
            {
                // SampleId,SigName,SigContrib,SigPercent
                final String[] items = line.split(DATA_DELIM, -1);
                String cancerType = items[0];
                String sigName = items[1];

                double[] percentileData = new double[PERCENTILE_COUNT];

                int startIndex = 2;

                for(int i = startIndex; i < items.length; ++i)
                {
                    double value = Double.parseDouble(items[i]);
                    percentileData[i - startIndex] = value;
                }

                Map<String,double[]> sigContribsMap = cancerSigContribsMap.get(cancerType);

                if(sigContribsMap == null)
                {
                    sigContribsMap = Maps.newHashMap();
                    cancerSigContribsMap.put(cancerType, sigContribsMap);
                }

                sigContribsMap.put(sigName, percentileData);
            }
        }
        catch (IOException e)
        {
            CUP_LOGGER.error("failed to read sig contrib percentile data file({}): {}", filename, e.toString());
        }
    }

}
