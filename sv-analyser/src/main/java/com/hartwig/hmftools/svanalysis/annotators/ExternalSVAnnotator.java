package com.hartwig.hmftools.svanalysis.annotators;

import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.svanalysis.types.SvVarData;
import com.google.common.collect.Maps;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;

public class ExternalSVAnnotator {

    private static final Logger LOGGER = LogManager.getLogger(ExternalSVAnnotator.class);

    private static int CSV_INDEX_SAMPLE_ID = 1;
    private static int CSV_INDEX_SV_ID = 2;
    private static int CSV_FIELD_COUNT = 9;
    private static int CSV_ID_FIELD_COUNT = 3;

    // annotations
    private static int FIELD_INDEX_PON_COUNT = 0;
    private static int FIELD_INDEX_PON_REGION_COUNT = 1;
    private static int FIELD_INDEX_LE_START = 2;
    private static int FIELD_INDEX_LE_END = 3;
    private static int FIELD_INDEX_FS_START = 4;
    private static int FIELD_INDEX_FS_END = 5;

    private Map<String, Map<String, ExternalSvData>> mSampleSvData;

    public ExternalSVAnnotator()
    {
        mSampleSvData = Maps.newHashMap();
    }

    public void loadFile(final String filename) {
        if (filename.isEmpty()) {
            return;
        }

        try {

            BufferedReader fileReader = new BufferedReader(new FileReader(filename));

            // read field names
            String line = fileReader.readLine();

            if (line == null) {
                LOGGER.error("Empty external SV annotations CSV file({})", filename);
                return;
            }

            int svCount = 0;

            while ((line = fileReader.readLine()) != null) {

                // parse CSV data
                String[] items = line.split(",");

                if (items.length != CSV_FIELD_COUNT) {
                    continue;
                }

                String sampleId = items[CSV_INDEX_SAMPLE_ID];

                if(!mSampleSvData.containsKey(sampleId))
                {
                    mSampleSvData.put(sampleId, new HashMap<String, ExternalSvData>());
                }

                Map<String, ExternalSvData> sampleDataMap = mSampleSvData.get(sampleId);

                final String svId = items[CSV_INDEX_SV_ID];
                ++svCount;

                List<String> dataValues = Lists.newArrayList();

                // now parse the required values
                dataValues.add(FIELD_INDEX_PON_COUNT, items[FIELD_INDEX_PON_COUNT+CSV_ID_FIELD_COUNT]);
                dataValues.add(FIELD_INDEX_PON_REGION_COUNT, items[FIELD_INDEX_PON_REGION_COUNT+CSV_ID_FIELD_COUNT]);
                dataValues.add(FIELD_INDEX_LE_START, items[FIELD_INDEX_LE_START+CSV_ID_FIELD_COUNT]);
                dataValues.add(FIELD_INDEX_LE_END, items[FIELD_INDEX_LE_END+CSV_ID_FIELD_COUNT]);
                dataValues.add(FIELD_INDEX_FS_START, items[FIELD_INDEX_FS_START+CSV_ID_FIELD_COUNT]);
                dataValues.add(FIELD_INDEX_FS_END, items[FIELD_INDEX_FS_END+CSV_ID_FIELD_COUNT]);

                sampleDataMap.put(svId, new ExternalSvData(svId, dataValues));
            }

            LOGGER.debug("loaded {} samples and {} SVs from external data file", mSampleSvData.size(), svCount);

        } catch (IOException exception) {
            LOGGER.error("Failed to read external SV annotations CSV file({})", filename);
        }
    }

    public boolean hasData() { return !mSampleSvData.isEmpty(); }

    public void setSVData(final String sampleId, SvVarData var) {

        if(!mSampleSvData.containsKey(sampleId))
            return;

        final Map<String, ExternalSvData> sampleDataMap = mSampleSvData.get(sampleId);

        final ExternalSvData svData = sampleDataMap.get(var.id());

        if(svData == null) {

            LOGGER.error("sv({}) external data not found", var.id());
            return;
        }

        final List<String> svValues = svData.getValues();

        var.setPonCount(Integer.parseInt(svValues.get(FIELD_INDEX_PON_COUNT)));
        var.setPonRegionCount(Integer.parseInt(svValues.get(FIELD_INDEX_PON_REGION_COUNT)));
        var.setLineElement(svValues.get(FIELD_INDEX_LE_START), true);
        var.setLineElement(svValues.get(FIELD_INDEX_LE_END), false);
        var.setFragileSites(svValues.get(FIELD_INDEX_FS_START), svValues.get(FIELD_INDEX_FS_END));
    }
}
