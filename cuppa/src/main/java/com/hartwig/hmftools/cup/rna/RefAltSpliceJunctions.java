package com.hartwig.hmftools.cup.rna;

import static com.hartwig.hmftools.common.utils.file.CommonFields.FLD_CHROMOSOME;
import static com.hartwig.hmftools.common.utils.file.CommonFields.FLD_GENE_ID;
import static com.hartwig.hmftools.common.utils.file.FileWriterUtils.createBufferedReader;
import static com.hartwig.hmftools.common.utils.MatrixFile.DEFAULT_MATRIX_DELIM;
import static com.hartwig.hmftools.common.utils.file.FileWriterUtils.closeBufferedWriter;
import static com.hartwig.hmftools.common.utils.file.FileWriterUtils.createBufferedWriter;
import static com.hartwig.hmftools.cup.CuppaConfig.CUP_LOGGER;
import static com.hartwig.hmftools.cup.CuppaConfig.DATA_DELIM;
import static com.hartwig.hmftools.cup.CuppaRefFiles.REF_FILE_ALT_SJ_CANCER;
import static com.hartwig.hmftools.common.cuppa.CategoryType.ALT_SJ;
import static com.hartwig.hmftools.common.cuppa.ClassifierType.ALT_SJ_COHORT;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.IOException;
import java.util.List;
import java.util.Map;

import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import com.hartwig.hmftools.common.utils.Matrix;
import com.hartwig.hmftools.common.cuppa.CategoryType;
import com.hartwig.hmftools.common.utils.config.ConfigBuilder;
import com.hartwig.hmftools.cup.common.NoiseRefCache;
import com.hartwig.hmftools.cup.common.SampleData;
import com.hartwig.hmftools.cup.common.SampleDataCache;
import com.hartwig.hmftools.cup.ref.RefClassifier;
import com.hartwig.hmftools.cup.ref.RefDataConfig;

import org.apache.commons.cli.CommandLine;

@Deprecated
public class RefAltSpliceJunctions implements RefClassifier
{
    private final RefDataConfig mConfig;
    private final SampleDataCache mSampleDataCache;

    public static final String FLD_POS_START = "PosStart";
    public static final String FLD_POS_END = "PosEnd";

    public RefAltSpliceJunctions(final RefDataConfig config, final SampleDataCache sampleDataCache)
    {
        mConfig = config;
        mSampleDataCache = sampleDataCache;
    }

    public CategoryType categoryType() { return ALT_SJ; }

    public static boolean requiresBuild(final RefDataConfig config)
    {
        return config.Categories.contains(ALT_SJ) || !config.AltSjMatrixFile.isEmpty();
    }

    @Override
    public boolean buildRefDataSets()
    {
        CUP_LOGGER.debug("loading RNA alternate splice-junction data");

        final List<String> cancerTypes = Lists.newArrayList();
        final Map<String,Integer> sampleIndexMap = Maps.newHashMap();
        final List<String> asjLocations = loadAltSjLocations(mConfig.AltSjMatrixFile);

        if(asjLocations == null)
            return false;

        int altSjSiteCount = asjLocations.size();
        final short[][] sampleFragCounts = loadSampleAltSjMatrixData(mConfig.AltSjMatrixFile, sampleIndexMap, altSjSiteCount);

        if(sampleFragCounts == null)
            return false;

        // alt-SJs in the columns, cancer types in the rows since the noise allocation expects this
        Matrix cancerFragCounts = new Matrix(mSampleDataCache.RefCancerSampleData.size(), altSjSiteCount);
        final double[][] cancerMatrixData = cancerFragCounts.getData();

        for(Map.Entry<String,List<SampleData>> entry : mSampleDataCache.RefCancerSampleData.entrySet())
        {
            final String refCancerType = entry.getKey();
            cancerTypes.add(refCancerType);
            int cancerIndex = cancerTypes.size() - 1;

            for(final SampleData sample : entry.getValue())
            {
                if(!sampleIndexMap.containsKey(sample.Id))
                    continue;

                Integer countsIndex = sampleIndexMap.get(sample.Id);

                if(countsIndex == null)
                {
                    CUP_LOGGER.warn("sample({}) missing alt-SJ data", sample.Id);
                    continue;
                }

                final short[] fragCounts = sampleFragCounts[countsIndex];

                for(int b = 0; b < fragCounts.length; ++b)
                {
                    cancerMatrixData[cancerIndex][b] += fragCounts[b];
                }
            }
        }

        final double[] altSjMedians = NoiseRefCache.generateMedianValues(cancerFragCounts);

        if(mConfig.NoiseAdjustments.hasNoiseAllocation(ALT_SJ_COHORT))
        {
            int noiseAllocation = mConfig.NoiseAdjustments.getNoiseAllocation(ALT_SJ_COHORT);
            CUP_LOGGER.debug("applying noise({}) to alt-SJ cohort counts", noiseAllocation);

            NoiseRefCache.applyNoise(cancerFragCounts, altSjMedians, noiseAllocation);
        }

        // no need to write medians unless to try out pairwise in the classifer

        CUP_LOGGER.debug("writing RNA alt-SJ cancer reference data");

        writeCancerAltSjMatrixData(cancerFragCounts, cancerTypes, asjLocations);
        return true;
    }

    public static final int ASJ_LOCATION_COL_COUNT = 4;

    private static List<String> loadAltSjLocations(final String filename)
    {
        final List<String> asjLocations = Lists.newArrayList();

        try
        {
            BufferedReader fileReader = createBufferedReader(filename);

            fileReader.readLine(); // skip header
            char delim = DEFAULT_MATRIX_DELIM.charAt(0);

            String line = null;
            while((line = fileReader.readLine()) != null)
            {
                // looks like: GeneId,Chromosome,PosStart,PosEnd, eg ENSG00000227232,1,14829,14930
                int charIndex = 0;
                int delims = 0;
                int lineLength = line.length();
                while(delims < 4 && charIndex < lineLength)
                {
                    if(line.charAt(charIndex) == delim)
                        ++delims;

                    ++charIndex;
                }

                if(delims != ASJ_LOCATION_COL_COUNT)
                {
                    CUP_LOGGER.error("invalid alt-SJ location header: {}", line);
                    return null;
                }

                asjLocations.add(line.substring(0, charIndex - 1));
             }
        }
        catch (IOException exception)
        {
            CUP_LOGGER.error("failed to read alt-SJ locations from matrix data file({}): {}", filename, exception.toString());
            return null;
        }

        return asjLocations;
    }

    public static short[][] loadSampleAltSjMatrixData(final String filename, final Map<String,Integer> sampleIndexMap, int altSjSiteCount)
    {
        // expect the matrix to start with columns GeneId,Chromosome,PosStart,PosEnd
        short[][] sampleCounts = null;

        try
        {
            BufferedReader fileReader = createBufferedReader(filename);

            String header = fileReader.readLine();
            final String[] columns = header.split(DEFAULT_MATRIX_DELIM, -1);

            for(int i = ASJ_LOCATION_COL_COUNT; i < columns.length; ++i)
            {
                sampleIndexMap.put(columns[i], i - ASJ_LOCATION_COL_COUNT);
            }

            int sampleCount = sampleIndexMap.size();
            sampleCounts = new short[sampleCount][altSjSiteCount];

            int altSjIndex = 0;
            int zeroCount = 0;
            int oneCount = 0;
            int maxShortCount = 0;

            String line = null;
            while((line = fileReader.readLine()) != null)
            {
                final String[] values = line.split(DATA_DELIM, -1);

                if(values.length != ASJ_LOCATION_COL_COUNT + sampleCount)
                {
                    CUP_LOGGER.error("invalid alt-SJ sample matrix column count({}) vs expected({} + {})",
                            values.length, ASJ_LOCATION_COL_COUNT, sampleCount);
                    return null;
                }

                int sampleIndex = 0;
                for(int i = ASJ_LOCATION_COL_COUNT; i < values.length; ++i)
                {
                    // zeros may be represented a an empty entry
                    int fragCountRaw = values[i].isEmpty() ? 0 : Integer.parseInt(values[i]);

                    short fragCount;
                    if(fragCountRaw > Short.MAX_VALUE)
                    {
                        ++maxShortCount;
                        fragCount = Short.MAX_VALUE;
                    }
                    else
                    {
                        fragCount = (short)fragCountRaw;
                    }

                    // data is transposed
                    sampleCounts[sampleIndex][altSjIndex] = fragCount;

                    if(fragCount == 0)
                        ++zeroCount;
                    else if(fragCount == 1)
                        ++oneCount;

                    ++sampleIndex;
                }

                ++altSjIndex;
            }

            CUP_LOGGER.info("loaded matrix(rows={} cols={}) from file({})", altSjSiteCount, sampleCount, filename);

            if(maxShortCount > 0)
            {
                CUP_LOGGER.warn("file({}) zeros({}) ones({}) capped {} short values", filename, zeroCount, oneCount, maxShortCount);
            }

        }
        catch (IOException exception)
        {
            CUP_LOGGER.error("failed to read matrix data file({}): {}", filename, exception.toString());
            return null;
        }

        return sampleCounts;
    }

    private void writeCancerAltSjMatrixData(
            final Matrix fragCountMatrix, final List<String> cancerTypes, final List<String> asjLocations)
    {
        try
        {
            final String filename = mConfig.OutputDir + REF_FILE_ALT_SJ_CANCER;
            BufferedWriter writer = createBufferedWriter(filename, false);

            writer.write(String.format("%s,%s,%s,%s", FLD_GENE_ID, FLD_CHROMOSOME, FLD_POS_START, FLD_POS_END));

            for(final String header : cancerTypes)
            {
                writer.write(String.format(",%s", header));
            }

            writer.newLine();

            final double[][] matrixData = fragCountMatrix.getData();

            for(int b = 0; b < fragCountMatrix.Cols; ++b) // columns are the alt-SJ locations
            {
                writer.write(String.format("%s", asjLocations.get(b)));

                for(int i = 0; i < fragCountMatrix.Rows; ++i)
                {
                    writer.write(String.format(",%.1f", matrixData[i][b]));
                }

                writer.newLine();
            }

            closeBufferedWriter(writer);
        }
        catch(IOException e)
        {
            CUP_LOGGER.error("failed to write ref RNA alt-SJ data: {}", e.toString());
        }
    }

}
