package com.hartwig.hmftools.cup.prep;

import static com.hartwig.hmftools.common.utils.file.FileDelimiters.TSV_DELIM;
import static com.hartwig.hmftools.cup.common.CupConstants.CUP_LOGGER;
import static com.hartwig.hmftools.cup.prep.DataItem.FLD_CATEGORY;
import static com.hartwig.hmftools.cup.prep.DataItem.FLD_KEY;
import static com.hartwig.hmftools.cup.prep.DataItem.FLD_SOURCE;
import static com.hartwig.hmftools.cup.prep.DataItem.FLD_VALUE;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;
import java.util.StringJoiner;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.utils.file.FileWriterUtils;

public class DataItemsIO
{
    private static final String[] INDEX_FIELDS = new String[] { FLD_SOURCE, FLD_CATEGORY, FLD_KEY };

    public static void writeDataItemList(List<DataItem> dataItems, String path)
    {
        try
        {
            CUP_LOGGER.info("Writing data to: " + path);

            BufferedWriter writer = FileWriterUtils.createBufferedWriter(path, false);

            List<String> headerFields = Lists.newArrayList(INDEX_FIELDS);
            headerFields.add(FLD_VALUE);

            String header = String.join(TSV_DELIM, headerFields);
            writer.write(header);
            writer.newLine();

            for(DataItem dataItem : dataItems)
            {
                String line = String.join(
                        TSV_DELIM,
                        dataItem.Index.Source.toString(),
                        dataItem.Index.Type.getAlias(),
                        dataItem.Index.Key,
                        dataItem.Value
                );

                writer.write(line);
                writer.newLine();
            }

            FileWriterUtils.closeBufferedWriter(writer);
        }
        catch(Exception e)
        {
            CUP_LOGGER.error("Failed to write features:");
            e.printStackTrace();
            System.exit(1);
        }
    }

    public static void writeDataItemMatrix(DataItemMatrix dataItemMatrix, String path, boolean append)
    {
        try
        {
            CUP_LOGGER.info("Writing data to: " + path);

            StringJoiner joiner = new StringJoiner(TSV_DELIM);
            BufferedWriter writer = FileWriterUtils.createBufferedWriter(path, append);

            if(!append)
            {
                for(String field : INDEX_FIELDS)
                    joiner.add(field);

                for(String sampleId : dataItemMatrix.SampleIds)
                    joiner.add(sampleId);

                String header = joiner.toString();
                writer.write(header);
                writer.newLine();
            }

            for(DataItem.Index index : dataItemMatrix.Indexes)
            {
                joiner = new StringJoiner(TSV_DELIM);

                joiner.add(index.Source.toString());
                joiner.add(index.Type.getAlias());
                joiner.add(index.Key);

                for(int sampleIndex = 0; sampleIndex < dataItemMatrix.nSamples(); sampleIndex++)
                {
                    String value = dataItemMatrix.get(index)[sampleIndex];
                    joiner.add(value);
                }

                writer.write(joiner.toString());
                writer.newLine();
            }

            FileWriterUtils.closeBufferedWriter(writer);
        }
        catch(Exception e)
        {
            CUP_LOGGER.error("Failed to write multi-sample feature matrix:");
            e.printStackTrace();
            System.exit(1);
        }
    }

    public static List<DataItem> readDataItemList(String path)
    {
        try
        {
            BufferedReader fileReader = FileWriterUtils.createBufferedReader(path);

            String[] header = fileReader.readLine().split(TSV_DELIM);
            int EXPECTED_HEADER_LENGTH = INDEX_FIELDS.length + 1;
            if(header.length > EXPECTED_HEADER_LENGTH)
                CUP_LOGGER.warn("Header has >{} fields. The input file might be a multi-sample file: {}", EXPECTED_HEADER_LENGTH, path);

            List<DataItem> dataItems = new ArrayList<>();
            String line;
            while((line = fileReader.readLine()) != null)
            {
                String[] rowValues = line.split(TSV_DELIM, -1);

                DataItem dataItem = new DataItem(
                        DataSource.valueOf(rowValues[0]),
                        ItemType.fromAlias(rowValues[1]),
                        rowValues[2],
                        rowValues[3]
                );

                dataItems.add(dataItem);
            }

            return dataItems;
        }
        catch(IOException e)
        {
            CUP_LOGGER.error("Failed to load data items from file: {}", path);
            e.printStackTrace();
            return null;
        }
    }

     public static DataItemMatrix readDataItemMatrix(String path)
     {
         try
         {
             BufferedReader fileReader = FileWriterUtils.createBufferedReader(path);

             String[] header = fileReader.readLine().split(TSV_DELIM);
             List<String> sampleIds = List.of(Arrays.copyOfRange(header, INDEX_FIELDS.length, header.length));

             Map<DataItem.Index, String[]> featureBySampleMatrix = new LinkedHashMap<>();
             String line;
             while((line = fileReader.readLine()) != null)
             {
                 String[] rowValues = line.split(TSV_DELIM, -1);

                 DataItem.Index index = new DataItem.Index(
                         DataSource.valueOf(rowValues[0]),
                         ItemType.fromAlias(rowValues[1]),
                         rowValues[2]
                 );

                 String[] featureValues = new String[sampleIds.size()];
                 int sampleIndex = 0;
                 for(int columnIndex = INDEX_FIELDS.length; columnIndex < header.length; columnIndex++)
                 {
                     String featureValue = rowValues[columnIndex];
                     featureValues[sampleIndex] = featureValue.equals("null") ? null : featureValue;
                     sampleIndex++;
                 }

                 featureBySampleMatrix.put(index, featureValues);
             }

             return new DataItemMatrix(sampleIds, featureBySampleMatrix);
         }
         catch(IOException e)
         {
             CUP_LOGGER.error("Failed to load data items from file: {}", path);
             e.printStackTrace();
             return null;
         }
     }
}
