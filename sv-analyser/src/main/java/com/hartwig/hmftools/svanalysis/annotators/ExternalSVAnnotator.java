package com.hartwig.hmftools.svanalysis.annotators;

import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.svanalysis.types.SvClusterData;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;

public class ExternalSVAnnotator {


    private List<String> mFieldNames;
    private String mFieldNamesStr;
    private Map<Integer, String> mIdValues;
    private String mEmptyValues;

    private static final Logger LOGGER = LogManager.getLogger(ExternalSVAnnotator.class);


    public ExternalSVAnnotator()
    {
        mFieldNames = Lists.newArrayList();
        mFieldNamesStr = "";
        mIdValues = new HashMap();
        mEmptyValues = "";
    }

    public void loadFile(final String filename)
    {
        if(filename.isEmpty())
            return;

        try {

            BufferedReader fileReader = new BufferedReader(new FileReader(filename));

            // read field names
            String line = fileReader.readLine();

            if(line == null)
            {
                LOGGER.error("Empty external SV annotations CSV file({})", filename);
                return;
            }

            String[] fieldNames = line.split(",");

            mFieldNames = Lists.newArrayList(fieldNames);

            // cache the field names, excluding the SV ID
            for(int i = 1; i < fieldNames.length; ++i) {
                mEmptyValues += ",";
                mFieldNamesStr += "," + fieldNames[i];
            }

            mFieldNamesStr = mFieldNamesStr.substring(1);

            while ((line = fileReader.readLine()) != null) {

                // parse CSV data
                String[] items = line.split(",");

                if(items.length != mFieldNames.size())
                    continue;

                int svId = Integer.parseInt(items[0]);

                String values = "";
                for(int i = 1; i < items.length; ++i)
                {
                    values += "," + items[i];
                }

                values = values.substring(1);

                mIdValues.put(svId, values);
            }

            LOGGER.debug("loaded {} external SV annotations", mIdValues.size());
        }
        catch(IOException exception)
        {
            LOGGER.error("Failed to read external SV annotations CSV file({})", filename);
        }
    }

    public boolean hasExternalData() { return !mFieldNames.isEmpty(); }
    public final String getFieldNames() { return mFieldNamesStr; }

    public final String getSVData(final SvClusterData var)
    {
        final String dataList = mIdValues.get(Integer.parseInt(var.id()));

        if(dataList == null)
            return mEmptyValues;

        return dataList;
    }
}
