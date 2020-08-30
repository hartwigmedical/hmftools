package com.hartwig.hmftools.common.sigs;

import static java.lang.Math.abs;

import static com.hartwig.hmftools.common.sigs.VectorUtils.sumVector;
import static com.hartwig.hmftools.common.utils.GenericDataCollection.GD_TYPE_DECIMAL;
import static com.hartwig.hmftools.common.utils.GenericDataCollection.GD_TYPE_INTEGER;
import static com.hartwig.hmftools.common.utils.GenericDataCollection.GD_TYPE_STRING;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.nio.file.Files;
import java.util.List;
import java.util.Map;

import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import com.hartwig.hmftools.common.utils.GenericDataCollection;
import com.hartwig.hmftools.common.utils.GenericDataLoader;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;

public class SigUtils
{
    public static final Logger SU_LOGGER = LogManager.getLogger(SigUtils.class);

    public static final String DEFAULT_MATRIX_DELIM = ",";

    public static void convertToPercentages(double[] counts)
    {
        double total = sumVector(counts);

        if(total <= 0)
            return;

        for(int i = 0; i < counts.length; ++i)
        {
            counts[i] /= total;
        }
    }

    public static final double[] calculateFittedCounts(final SigMatrix signatures, final double[] allocations)
    {
        double[] fittedCounts = new double[signatures.Rows];

        for(int transId = 0; transId < signatures.Cols; ++transId)
        {
            double allocation = allocations[transId];

            for(int catId = 0; catId < signatures.Rows; ++catId)
            {
                fittedCounts[catId] += allocation * signatures.get(catId, transId);
            }
        }

        return fittedCounts;
    }

    public static SigResiduals calcResiduals(final double[] counts, final double[] fittedCounts, double totalCounts)
    {
        SigResiduals residuals = new SigResiduals();

        for(int catId = 0; catId < counts.length; ++catId)
        {
            double diff = fittedCounts[catId] - counts[catId];
            residuals.Total += abs(diff);

            if(diff > 0)
                residuals.Excess += diff;
        }

        residuals.Percent = residuals.Total / totalCounts;
        return residuals;
    }

    public static double calcLinearLeastSquares(final double[] params, final double[] data)
    {
        if(data.length != params.length)
            return 0;

        // returns the best ratio applying the params to the data assuming direct ratio (ie line through origin)
        double paramTotal = 0;
        double multTotal = 0;

        for(int i = 0; i < data.length; ++i)
        {
            paramTotal += params[i] * params[i];
            multTotal += params[i] * data[i];
        }

        return paramTotal > 0 ? multTotal/paramTotal : 0;
    }

    public static void copyMatrix(final double[][] source, final double[][] dest)
    {
        if(source.length != dest.length)
            return;

        int rows = source.length;
        int cols = source[0].length;

        for(int i = 0; i < rows; ++i)
        {
            for (int j = 0; j < cols; ++j)
            {
                dest[i][j] = source[i][j];
            }
        }
    }

    public static double calcAbsDiffs(final double[] set1, final double[] set2)
    {
        if(set1.length != set2.length)
            return 0;

        double diffTotal = 0;
        for(int i = 0; i < set1.length; ++i)
        {
            diffTotal += abs(set1[i] - set2[i]);
        }

        return diffTotal;
    }

    public static SigMatrix createMatrixFromListData(final List<List<Double>> dataSet)
    {
        if(dataSet.isEmpty())
            return null;

        int rows = dataSet.size();
        int cols = dataSet.get(0).size();

        // create and populate the matrix
        SigMatrix matrix = new SigMatrix(rows, cols);

        double[][] matrixData = matrix.getData();

        for(int i = 0; i < rows; ++i)
        {
            for (int j = 0; j < cols; ++j)
            {
                matrixData[i][j] = dataSet.get(i).get(j);
            }
        }

        return matrix;
    }

    public static SigMatrix loadMatrixDataFile(final String filename, final Map<String,Integer> columnDataIndex, final List<String> ignoreFields)
    {
        final List<String> columnNames = Lists.newArrayList();

        SigMatrix sigMatrix = loadMatrixDataFile(filename, columnNames, ignoreFields);

        for(int c = 0; c < columnNames.size(); ++c)
        {
            columnDataIndex.put(columnNames.get(c), c);
        }

        return sigMatrix;
    }

    public static SigMatrix loadMatrixDataFile(final String filename, final List<String> columnNames)
    {
        return loadMatrixDataFile(filename, columnNames, null);
    }

    public static SigMatrix loadMatrixDataFile(final String filename, final List<String> columnNames, final List<String> ignoreFields)
    {
        try
        {
            final List<String> fileData = Files.readAllLines(new File(filename).toPath());

            // read field names
            if(fileData.size() <= 1)
            {
                SU_LOGGER.error("empty data CSV file({})", filename);
                return null;
            }

            final String header = fileData.get(0);
            columnNames.addAll(Lists.newArrayList(header.split(DEFAULT_MATRIX_DELIM, -1)));
            fileData.remove(0);

            List<Integer> ignoreCols = Lists.newArrayList();

            if(ignoreFields != null)
            {
                for(int i = 0; i < columnNames.size(); ++i)
                {
                    if(ignoreFields.contains(columnNames.get(i)))
                    {
                        ignoreCols.add(i);

                        if(ignoreCols.size() == ignoreFields.size())
                            break;
                    }
                }

                ignoreFields.forEach(x -> columnNames.remove(x));
            }

            int colCount = columnNames.size();
            int rowCount = fileData.size();
            int invalidRowCount = 0;

            SigMatrix matrix = new SigMatrix(rowCount, colCount);
            final double[][] matrixData = matrix.getData();

            for(int r = 0; r < rowCount; ++r)
            {
                final String line = fileData.get(r);
                final String[] items = line.split(DEFAULT_MATRIX_DELIM, -1);

                if(items.length != colCount + ignoreCols.size())
                {
                    ++invalidRowCount;
                    continue;
                }

                int c = 0;
                for(int i = 0; i < items.length; ++i)
                {
                    if(ignoreCols.contains(i))
                        continue;

                    matrixData[r][c] = Double.parseDouble(items[i]);
                    ++c;
                }
            }

            SU_LOGGER.info("loaded matrix(rows={} cols={}) from file({}), invalid rows({})", rowCount, colCount, filename, invalidRowCount);

            return matrix;

        }
        catch (IOException exception)
        {
            SU_LOGGER.error("failed to read matrix data file({}): {}", filename, exception.toString());
            return null;
        }
    }

    public static void writeMatrixData(
            final BufferedWriter writer, final List<String> headers, final SigMatrix matrix, boolean asInt) throws IOException
    {
        if(headers != null)
        {
            int i = 0;
            for (; i < headers.size() - 1; ++i)
            {
                writer.write(String.format("%s,", headers.get(i)));
            }
            writer.write(String.format("%s", headers.get(i)));

            writer.newLine();
        }

        final double[][] sigData = matrix.getData();

        for(int i = 0; i < matrix.Rows; ++i)
        {
            for(int j = 0; j < matrix.Cols; ++j)
            {
                if(asInt)
                    writer.write(String.format("%.0f", sigData[i][j]));
                else
                    writer.write(String.format("%.6f", sigData[i][j]));

                if(j < matrix.Cols-1)
                    writer.write(String.format(",", sigData[i][j]));
            }

            writer.newLine();
        }

    }

    public static void writeMatrixData(final BufferedWriter writer, final SigMatrix matrix, boolean asInt) throws IOException
    {
        writeMatrixData(writer, null, matrix, asInt);
    }


}
