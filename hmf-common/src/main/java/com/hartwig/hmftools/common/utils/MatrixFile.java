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

public final class MatrixFile
{
    public static final String DEFAULT_MATRIX_DELIM = ",";

    private static final Logger LOGGER = LogManager.getLogger(MatrixFile.class);

    public static Matrix loadMatrixDataFile(
            final String filename, final Map<String,Integer> columnDataIndex, final List<String> ignoreFields, boolean transpose)
    {
        final List<String> columnNames = Lists.newArrayList();

        Matrix sigMatrix = loadMatrixDataFile(filename, columnNames, ignoreFields, transpose);

        for(int c = 0; c < columnNames.size(); ++c)
        {
            columnDataIndex.put(columnNames.get(c), c);
        }

        return sigMatrix;
    }

    public static Matrix loadMatrixDataFile(final String filename, final List<String> columnNames, boolean transpose)
    {
        return loadMatrixDataFile(filename, columnNames, null, transpose);
    }

    public static Matrix loadMatrixDataFile(
            final String filename, final List<String> columnNames, final List<String> ignoreFields, boolean transpose)
    {
        try
        {
            final List<String> fileData = Files.readAllLines(new File(filename).toPath());

            Matrix matrix = loadMatrixDataFile(fileData, columnNames, ignoreFields, transpose);

            LOGGER.info("loaded matrix(rows={} cols={}) from file({})", matrix.Rows, matrix.Cols, filename);

            return matrix;
        }
        catch (IOException exception)
        {
            LOGGER.error("failed to read matrix data file({}): {}", filename, exception.toString());
            return null;
        }
    }

    public static Matrix loadMatrixDataFile(
            final List<String> fileData, final List<String> columnNames, final List<String> ignoreFields, boolean transpose)
    {
        // read field names
        if(fileData.size() <= 1)
        {
            LOGGER.error("empty matrix data");
            return null;
        }

        final String header = fileData.get(0);
        fileData.remove(0);
        columnNames.addAll(Lists.newArrayList(header.split(DEFAULT_MATRIX_DELIM, -1)));

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

        Matrix matrix = transpose ? new Matrix(colCount, rowCount) : new Matrix(rowCount, colCount);
        final double[][] matrixData = matrix.getData();

        for(int rowIndex = 0; rowIndex < rowCount; ++rowIndex)
        {
            final String line = fileData.get(rowIndex);
            final String[] items = line.split(DEFAULT_MATRIX_DELIM, -1);

            if(items.length != colCount + ignoreCols.size())
            {
                ++invalidRowCount;
                continue;
            }

            int colIndex = 0;
            for(int i = 0; i < items.length; ++i)
            {
                if(ignoreCols.contains(i))
                    continue;

                if(transpose)
                    matrixData[colIndex][rowIndex] = Double.parseDouble(items[i]);
                else
                    matrixData[rowIndex][colIndex] = Double.parseDouble(items[i]);

                ++colIndex;
            }
        }

        if(invalidRowCount > 0)
        {
            LOGGER.warn("matrix(rows={} cols={}) has {} invalid rows", matrix.Rows, matrix.Cols, invalidRowCount);
        }

        return matrix;
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
