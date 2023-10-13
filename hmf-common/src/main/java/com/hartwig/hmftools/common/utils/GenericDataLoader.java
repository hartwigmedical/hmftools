package com.hartwig.hmftools.common.utils;

import static com.hartwig.hmftools.common.utils.GenericDataCollection.GD_TYPE_DECIMAL;
import static com.hartwig.hmftools.common.utils.GenericDataCollection.GD_TYPE_INTEGER;
import static com.hartwig.hmftools.common.utils.GenericDataCollection.GD_TYPE_STRING;
import static com.hartwig.hmftools.common.utils.GenericDataCollection.isValidType;

import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;
import java.util.List;

import com.google.common.collect.Lists;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;

public final class GenericDataLoader
{
    private static final Logger LOGGER = LogManager.getLogger(GenericDataLoader.class);

    public static final String DEFAULT_DELIM = ",";

    public static GenericDataCollection loadFile(final String filename)
    {
        return loadFile(filename, GD_TYPE_DECIMAL, DEFAULT_DELIM, Lists.newArrayList());
    }

    public static GenericDataCollection loadFile(final String filename, int dataType, final String delimiter,
            final List<String> ignoreFields)
    {
        if(filename == null || filename.isEmpty())
        {
            return null;
        }

        if(!isValidType(dataType))
        {
            return null;
        }

        try
        {

            BufferedReader fileReader = new BufferedReader(new FileReader(filename));

            // read field names
            String line = fileReader.readLine();

            if(line == null)
            {
                LOGGER.error("Empty data CSV file({})", filename);
                return null;
            }

            int invalidRowCount = 0;

            GenericDataCollection collection = new GenericDataCollection(dataType);

            final String[] fieldNames = line.split(delimiter);
            int reqFieldCount = fieldNames.length;
            final List<String> headerNames = Lists.newArrayList(fieldNames);

            List<Integer> ignoreCols = Lists.newArrayList();
            for(int i = 0; i < headerNames.size(); ++i)
            {
                if(ignoreFields.contains(headerNames.get(i)))
                {
                    ignoreCols.add(i);
                }
            }

            ignoreFields.forEach(x -> headerNames.remove(x));

            collection.setFieldNames(headerNames);

            while((line = fileReader.readLine()) != null)
            {
                // parse CSV data
                String[] items = dataType == GD_TYPE_STRING ? line.split(",", -1) : line.split(","); // to include empty values

                if(items.length != reqFieldCount)
                {
                    ++invalidRowCount;
                    continue;
                }

                if(dataType == GD_TYPE_DECIMAL)
                {
                    List<Double> dataValues = Lists.newArrayList();

                    for(int i = 0; i < items.length; ++i)
                    {
                        if(ignoreCols.contains(i))
                        {
                            continue;
                        }

                        dataValues.add(Double.parseDouble(items[i]));
                    }

                    collection.addDecimalValues(dataValues);
                }
                else if(dataType == GD_TYPE_INTEGER)
                {
                    List<Integer> dataValues = Lists.newArrayList();

                    for(int i = 0; i < items.length; ++i)
                    {
                        if(ignoreCols.contains(i))
                        {
                            continue;
                        }

                        dataValues.add(Integer.parseInt(items[i]));
                    }

                    collection.addIntValues(dataValues);
                }
                else if(dataType == GD_TYPE_STRING)
                {
                    List<String> dataValues = Lists.newArrayList();

                    for(int i = 0; i < items.length; ++i)
                    {
                        if(ignoreCols.contains(i))
                        {
                            continue;
                        }

                        dataValues.add(items[i]);
                    }

                    collection.addStringValues(dataValues);
                }
            }

            LOGGER.debug("loaded {} data sets, invalid rows({})", collection.getDataCount(), invalidRowCount);

            return collection;

        }
        catch(IOException exception)
        {
            LOGGER.error("failed to read data file({}): {}", filename, exception.toString());
            return null;
        }
    }

}
