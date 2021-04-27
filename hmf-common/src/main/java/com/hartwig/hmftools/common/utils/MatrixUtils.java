package com.hartwig.hmftools.common.utils;

import java.io.BufferedWriter;
import java.io.File;
import java.io.IOException;
import java.nio.file.Files;
import java.util.List;
import java.util.Map;

import com.google.common.collect.Lists;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;

public class MatrixUtils
{
    private static final Logger LOGGER = LogManager.getLogger(MatrixUtils.class);

    public static final String DEFAULT_MATRIX_DELIM = ",";

    public static Matrix createMatrixFromListData(final List<List<Double>> dataSet)
    {
        if(dataSet.isEmpty())
            return null;

        int rows = dataSet.size();
        int cols = dataSet.get(0).size();

        // create and populate the matrix
        Matrix matrix = new Matrix(rows, cols);

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

    public static Matrix loadMatrixDataFile(final String filename, final Map<String,Integer> columnDataIndex, final List<String> ignoreFields)
    {
        final List<String> columnNames = Lists.newArrayList();

        Matrix sigMatrix = loadMatrixDataFile(filename, columnNames, ignoreFields);

        for(int c = 0; c < columnNames.size(); ++c)
        {
            columnDataIndex.put(columnNames.get(c), c);
        }

        return sigMatrix;
    }

    public static Matrix loadMatrixDataFile(final String filename, final List<String> columnNames)
    {
        return loadMatrixDataFile(filename, columnNames, null);
    }

    public static Matrix loadMatrixDataFile(final String filename, final List<String> columnNames, final List<String> ignoreFields)
    {
        try
        {
            final List<String> fileData = Files.readAllLines(new File(filename).toPath());

            // read field names
            if(fileData.size() <= 1)
            {
                LOGGER.error("empty data CSV file({})", filename);
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

            Matrix matrix = new Matrix(rowCount, colCount);
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

            LOGGER.info("loaded matrix(rows={} cols={}) from file({}) {}",
                    rowCount, colCount, filename, invalidRowCount > 0 ? String.format(", invalidRowCount(%d)", invalidRowCount) : "");

            return matrix;

        }
        catch (IOException exception)
        {
            LOGGER.error("failed to read matrix data file({}): {}", filename, exception.toString());
            return null;
        }
    }

    public static void writeMatrixData(
            final BufferedWriter writer, final List<String> headers, final Matrix matrix, boolean asInt) throws IOException
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

    public static void writeMatrixData(final BufferedWriter writer, final Matrix matrix, boolean asInt) throws IOException
    {
        writeMatrixData(writer, null, matrix, asInt);
    }

}
