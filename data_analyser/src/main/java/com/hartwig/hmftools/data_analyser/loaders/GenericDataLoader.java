package com.hartwig.hmftools.data_analyser.loaders;

import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;
import java.util.List;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.data_analyser.types.GenericDataCollection;

import org.apache.logging.log4j.Level;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;

public class GenericDataLoader {

    private static final Logger LOGGER = LogManager.getLogger(GenericDataLoader.class);

    public GenericDataLoader()
    {
    }

    public static GenericDataCollection loadFile(final String filename) {

        if (filename == null || filename.isEmpty()) {
            return null;
        }

        try {

            BufferedReader fileReader = new BufferedReader(new FileReader(filename));

            // read field names
            String line = fileReader.readLine();

            if (line == null) {
                LOGGER.error("Empty data CSV file({})", filename);
                return null;
            }

            GenericDataCollection collection = new GenericDataCollection();

            String[] fieldNames = line.split(",");
            List<String> strList = Lists.newArrayList(fieldNames);

            collection.setFieldNames(strList);

            while ((line = fileReader.readLine()) != null) {

                // parse CSV data
                String[] items = line.split(",");

                if(items.length == 0)
                {
                    continue;
                }

                List<Double> dataValues = Lists.newArrayList();

                for (int i = 0; i < items.length; ++i) {
                    dataValues.add(Double.parseDouble(items[i]));
                }

                collection.addDataValues(dataValues);
            }

            LOGGER.debug("loaded {} data sets", collection.getData().size());

            return collection;

        } catch (IOException exception) {

            LOGGER.error("failed to read data file({})", filename);
            return null;
        }
    }

    // public boolean hasData() { return !mSampleSvData.isEmpty(); }

}
