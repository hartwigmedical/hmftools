package com.hartwig.hmftools.cup.rna;

import static com.hartwig.hmftools.common.rna.RnaCommon.FLD_CHROMOSOME;
import static com.hartwig.hmftools.common.rna.RnaCommon.FLD_GENE_ID;
import static com.hartwig.hmftools.common.utils.FileWriterUtils.createBufferedReader;
import static com.hartwig.hmftools.common.utils.MatrixUtils.DEFAULT_MATRIX_DELIM;
import static com.hartwig.hmftools.common.utils.FileWriterUtils.closeBufferedWriter;
import static com.hartwig.hmftools.common.utils.FileWriterUtils.createBufferedWriter;
import static com.hartwig.hmftools.cup.CuppaConfig.CUP_LOGGER;
import static com.hartwig.hmftools.cup.CuppaConfig.DATA_DELIM;
import static com.hartwig.hmftools.cup.CuppaRefFiles.REF_FILE_ALT_SJ_CANCER;
import static com.hartwig.hmftools.cup.common.CategoryType.ALT_SJ;
import static com.hartwig.hmftools.cup.common.CupCalcs.addPanCancerNoise;
import static com.hartwig.hmftools.cup.common.CupConstants.ALT_SJ_NOISE_ALLOCATION;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.IOException;
import java.nio.file.Files;
import java.util.List;
import java.util.Map;
import java.util.StringJoiner;

import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import com.hartwig.hmftools.common.utils.Matrix;
import com.hartwig.hmftools.cup.common.CategoryType;
import com.hartwig.hmftools.cup.common.SampleData;
import com.hartwig.hmftools.cup.common.SampleDataCache;
import com.hartwig.hmftools.cup.ref.RefClassifier;
import com.hartwig.hmftools.cup.ref.RefDataConfig;

import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.Options;
import org.jetbrains.annotations.NotNull;

public class RefAltSpliceJunctions implements RefClassifier
{
    private final RefDataConfig mConfig;
    private final SampleDataCache mSampleDataCache;

    private final int mNoiseAllocation;

    public static final String FLD_POS_START = "PosStart";
    public static final String FLD_POS_END = "PosEnd";

    public static final String ALT_SJ_NOISE_ALLOC = "alt_sj_noise_alloc";

    public RefAltSpliceJunctions(final RefDataConfig config, final SampleDataCache sampleDataCache, final CommandLine cmd)
    {
        mConfig = config;
        mSampleDataCache = sampleDataCache;

        mNoiseAllocation = cmd != null ? Integer.parseInt(cmd.getOptionValue(ALT_SJ_NOISE_ALLOC, String.valueOf(ALT_SJ_NOISE_ALLOCATION))) : 0;
    }

    public static void addCmdLineArgs(@NotNull Options options)
    {
        options.addOption(ALT_SJ_NOISE_ALLOC, true, "Alt-SJ noise allocation");
    }

    public CategoryType categoryType() { return ALT_SJ; }

    public static boolean requiresBuild(final RefDataConfig config)
    {
        return config.Categories.contains(ALT_SJ) || !config.AltSjMatrixFile.isEmpty();
    }

    public void buildRefDataSets()
    {
        CUP_LOGGER.debug("loading RNA alternate splice-junction data");

        final List<String> cancerTypes = Lists.newArrayList();
        final Map<String,Integer> sampleIndexMap = Maps.newHashMap();
        final List<String> asjLocations = loadAltSjLocations(mConfig.AltSjMatrixFile);

        if(asjLocations == null)
            return;

        int altSjSiteCount = asjLocations.size();
        final short[][] sampleFragCounts = loadSampleAltSjMatrixData(mConfig.AltSjMatrixFile, sampleIndexMap, altSjSiteCount);

        if(sampleFragCounts == null)
            return;

        // columns contain the alt-SJ locations
        Matrix cancerFragCounts = new Matrix(altSjSiteCount, mSampleDataCache.RefCancerSampleData.size());
        final double[][] cancerMatrixData = cancerFragCounts.getData();

        for(Map.Entry<String,List<SampleData>> entry : mSampleDataCache.RefCancerSampleData.entrySet())
        {
            final String refCancerType = entry.getKey();
            cancerTypes.add(refCancerType);
            int cancerIndex = cancerTypes.size() - 1;

            for(final SampleData sample : entry.getValue())
            {
                Integer countsIndex = sampleIndexMap.get(sample.Id);

                if(countsIndex == null)
                {
                    CUP_LOGGER.warn("sample({}) missing alt-SJ data", sample.Id);
                    continue;
                }

                final short[] fragCounts = sampleFragCounts[countsIndex];

                for(int b = 0; b < fragCounts.length; ++b)
                {
                    cancerMatrixData[b][cancerIndex] += fragCounts[b];
                }
            }
        }

        if(mNoiseAllocation > 0)
        {
            CUP_LOGGER.debug("applying alt-SJ noise({}) to cancer matrix", mNoiseAllocation);
            addPanCancerNoise(cancerFragCounts, mNoiseAllocation);
        }

        CUP_LOGGER.debug("writing RNA alt-SJ cancer reference data");

        writeMatrixData(cancerFragCounts, cancerTypes, asjLocations);
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

    // may be used again if the number of sites shrinks considerably
    public static Matrix loadSampleAltSjMatrixData(
            final String filename, final Map<String,Integer> sampleIndexMap, final List<String> asjLocations)
    {
        // expect the matrix to start with columns GeneId,Chromosome,PosStart,PosEnd
        Matrix sampleMatrix = null;

        try
        {
            final List<String> fileData = Files.readAllLines(new File(filename).toPath());

            // read field names
            if(fileData.size() <= 1)
            {
                CUP_LOGGER.error("empty data CSV file({})", filename);
                return sampleMatrix;
            }

            final String header = fileData.get(0);
            final String[] columns = header.split(DEFAULT_MATRIX_DELIM, -1);
            fileData.remove(0);

            for(int i = ASJ_LOCATION_COL_COUNT; i < columns.length; ++i)
            {
                sampleIndexMap.put(columns[i], i - ASJ_LOCATION_COL_COUNT);
            }

            int sampleCount = sampleIndexMap.size();
            int altSjCount = fileData.size();

            sampleMatrix = new Matrix(sampleCount, altSjCount);
            final double[][] matrixData = sampleMatrix.getData();

            for(int r = 0; r < altSjCount; ++r)
            {
                final String line = fileData.get(r);
                final String[] values = line.split(DEFAULT_MATRIX_DELIM, -1);

                if(values.length != ASJ_LOCATION_COL_COUNT + sampleCount)
                {
                    CUP_LOGGER.error("invalid alt-SJ sample matrix column count({}) vs expeted({} + {})",
                            values.length, ASJ_LOCATION_COL_COUNT, sampleCount);
                    return null;
                }

                StringJoiner asjLocation = new StringJoiner(DEFAULT_MATRIX_DELIM);
                for(int i = 0; i < ASJ_LOCATION_COL_COUNT; ++i)
                {
                    asjLocation.add(values[i]);
                }

                asjLocations.add(asjLocation.toString());

                int c = 0;
                for(int i = ASJ_LOCATION_COL_COUNT; i < values.length; ++i)
                {
                    // data transposed
                    int fragCount = Integer.parseInt(values[i]);
                    matrixData[c][r] = fragCount;
                    ++c;
                }
            }

            CUP_LOGGER.info("loaded matrix(rows={} cols={}) from file({})", altSjCount, sampleCount, filename);
        }
        catch (IOException exception)
        {
            CUP_LOGGER.error("failed to read matrix data file({}): {}", filename, exception.toString());
            return null;
        }

        return sampleMatrix;
    }

    private void writeMatrixData(
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

            for(int i = 0; i < fragCountMatrix.Rows; ++i) // rows are the alt-SJ locations
            {
                writer.write(String.format("%s", asjLocations.get(i)));

                for(int j = 0; j < fragCountMatrix.Cols; ++j)
                {
                    writer.write(String.format(",%.1f", matrixData[i][j]));
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
