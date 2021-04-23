package com.hartwig.hmftools.cup.rna;

import static com.hartwig.hmftools.common.rna.RnaCommon.FLD_CHROMOSOME;
import static com.hartwig.hmftools.common.rna.RnaCommon.FLD_GENE_ID;
import static com.hartwig.hmftools.common.utils.MatrixUtils.DEFAULT_MATRIX_DELIM;
import static com.hartwig.hmftools.common.utils.FileWriterUtils.closeBufferedWriter;
import static com.hartwig.hmftools.common.utils.FileWriterUtils.createBufferedWriter;
import static com.hartwig.hmftools.cup.CuppaConfig.CUP_LOGGER;
import static com.hartwig.hmftools.cup.CuppaRefFiles.REF_FILE_ALT_SJ_CANCER;
import static com.hartwig.hmftools.cup.common.CategoryType.ALT_SJ;

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

import org.apache.commons.cli.Options;
import org.jetbrains.annotations.NotNull;

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

    public static void addCmdLineArgs(@NotNull Options options)
    {
        // options.addOption(TPM_LOG_CUTOFF, true, "RNA TPM cut-off in log scale (default=0, not applied)");
    }

    public CategoryType categoryType() { return ALT_SJ; }

    public static boolean requiresBuild(final RefDataConfig config) { return !config.RefAltSjFile.isEmpty(); }

    public void buildRefDataSets()
    {
        CUP_LOGGER.debug("loading RNA gene expression data");

        final List<String> cancerTypes = Lists.newArrayList();
        final Map<String,Integer> sampleIndexMap = Maps.newHashMap();
        final List<String> asjLocations = Lists.newArrayList();
        final Matrix sampleFragCounts = loadSampleAltSjMatrixData(mConfig.RefAltSjFile, sampleIndexMap, asjLocations);

        if(sampleFragCounts == null)
            return;

        // columns contain the alt-SJ locations
        Matrix cancerFragCounts = new Matrix(sampleFragCounts.Cols, mSampleDataCache.RefCancerSampleData.size());
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

                final double[] fragCounts = sampleFragCounts.getRow(countsIndex);

                for(int b = 0; b < fragCounts.length; ++b)
                {
                    cancerMatrixData[b][cancerIndex] += fragCounts[b];
                }
            }
        }

        CUP_LOGGER.debug("writing RNA alt-SJ cancer reference data");

        writeMatrixData(cancerFragCounts, cancerTypes, asjLocations);
    }

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

            int asjLocationColCount = 4;

            for(int i = asjLocationColCount; i < columns.length; ++i)
            {
                sampleIndexMap.put(columns[i], i - asjLocationColCount);
            }

            int sampleCount = sampleIndexMap.size();
            int altSjCount = fileData.size();
            int invalidRowCount = 0;

            sampleMatrix = new Matrix(sampleCount, altSjCount);
            final double[][] matrixData = sampleMatrix.getData();

            for(int r = 0; r < altSjCount; ++r)
            {
                final String line = fileData.get(r);
                final String[] items = line.split(DEFAULT_MATRIX_DELIM, -1);

                if(items.length != asjLocationColCount + sampleCount)
                {
                    ++invalidRowCount;
                    continue;
                }

                StringJoiner asjLocation = new StringJoiner(DEFAULT_MATRIX_DELIM);
                for(int i = 0; i < asjLocationColCount; ++i)
                {
                    asjLocation.add(items[i]);
                }

                asjLocations.add(asjLocation.toString());

                int c = 0;
                for(int i = asjLocationColCount; i < items.length; ++i)
                {
                    // data transposed
                    matrixData[c][r] = Double.parseDouble(items[i]);
                    ++c;
                }
            }

            CUP_LOGGER.info("loaded matrix(rows={} cols={}) from file({}) {}",
                    altSjCount, sampleCount, filename, invalidRowCount > 0 ? String.format(", invalidRowCount(%d)", invalidRowCount) : "");
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
                    writer.write(String.format(",%.0f", matrixData[i][j]));
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
