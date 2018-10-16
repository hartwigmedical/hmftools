package com.hartwig.hmftools.svanalysis.annotators;

import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import com.hartwig.hmftools.common.variant.structural.StructuralVariantType;
import com.hartwig.hmftools.svanalysis.types.SvVarData;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;

public class ExtDataLinker {

    private static final Logger LOGGER = LogManager.getLogger(ExtDataLinker.class);

    private static int CSV_INDEX_SAMPLE_ID = 0;
    private static int CSV_INDEX_SV_ID = 1;
    private static int CSV_ID_FIELD_COUNT = 2;

    // SV positional fields
    private static int FIELD_INDEX_TYPE = 0;
    private static int FIELD_INDEX_CHR_START = 1;
    private static int FIELD_INDEX_POS_START = 2;
    private static int FIELD_INDEX_ORIENT_START = 3;
    private static int FIELD_INDEX_CHR_END = 4;
    private static int FIELD_INDEX_POS_END = 5;
    private static int FIELD_INDEX_ORIENT_END = 6;

    private static int CSV_MIN_FIELD_COUNT = 9;
    private static int FIELD_INDEX_DATA = 0;

    private Map<String, Map<String, ExternalSvData>> mSampleSvData;
    private boolean mWriteOutput;
    private boolean mAnnotateSV;

    public ExtDataLinker()
    {
        mSampleSvData = Maps.newHashMap();
        mWriteOutput = true;
        mAnnotateSV = false;
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

                if (items.length < CSV_MIN_FIELD_COUNT) {
                    continue;
                }

                String sampleId = items[CSV_INDEX_SAMPLE_ID];

                if(!mSampleSvData.containsKey(sampleId))
                {
                    mSampleSvData.put(sampleId, new HashMap<String, ExternalSvData>());
                }

                Map<String, ExternalSvData> sampleDataMap = mSampleSvData.get(sampleId);

                // final String svId = items[CSV_INDEX_SV_ID];

                String svPosId = items[CSV_ID_FIELD_COUNT+FIELD_INDEX_TYPE];
                svPosId +=  "_" + items[CSV_ID_FIELD_COUNT+FIELD_INDEX_CHR_START];
                svPosId +=  "_" + items[CSV_ID_FIELD_COUNT+FIELD_INDEX_POS_START];
                svPosId +=  "_" + items[CSV_ID_FIELD_COUNT+FIELD_INDEX_ORIENT_START];
                svPosId +=  "_" + items[CSV_ID_FIELD_COUNT+FIELD_INDEX_CHR_END];
                svPosId +=  "_" + items[CSV_ID_FIELD_COUNT+FIELD_INDEX_POS_END];
                svPosId +=  "_" + items[CSV_ID_FIELD_COUNT+FIELD_INDEX_ORIENT_END];

                ++svCount;

                List<String> dataValues = Lists.newArrayList();

                for(int i = 0; i < items.length - CSV_MIN_FIELD_COUNT; ++i) {
                    // now parse the required values
                    dataValues.add(items[CSV_MIN_FIELD_COUNT + i]);
                }

                sampleDataMap.put(svPosId, new ExternalSvData(svPosId, dataValues));
            }

            LOGGER.debug("loaded {} samples and {} SVs from external data file", mSampleSvData.size(), svCount);

        } catch (IOException exception) {

            LOGGER.error("Failed to read external SV annotations CSV file({})", filename);
        }
    }

    public boolean hasData() { return !mSampleSvData.isEmpty(); }

    public void setSVData(final String sampleId, List<SvVarData> svVarData) {

        if(!mSampleSvData.containsKey(sampleId))
            return;

        final Map<String, ExternalSvData> sampleDataMap = mSampleSvData.get(sampleId);

        for(SvVarData var : svVarData) {

            String svPosId = var.type().toString();
            svPosId += "_" + var.chromosome(true);
            svPosId += "_" + var.position(true);
            svPosId += "_" + var.orientation(true);
            svPosId += "_" + var.chromosome(false);
            svPosId += "_" + var.position(false);
            svPosId += "_" + var.orientation(false);

            final ExternalSvData svData = sampleDataMap.get(svPosId);

            if(var.type() == StructuralVariantType.INS)
            {
                LOGGER.debug("sv({} posId={}) INS found", var.id(), svPosId);
            }

            if (svData == null) {

                // not every SV will have an entry - depends on the type of data
                LOGGER.debug("sv({} posId={}) external data not found", var.id(), svPosId);
                continue;
            }

            final List<String> svValues = svData.getValues();

            if(svValues.isEmpty())
                continue;

            String valuesStr = svValues.get(0);

            for(int i = 1; i < svValues.size(); ++i)
                valuesStr += "," + svValues.get(i);

            if(mWriteOutput)
            {
                LOGGER.info("EXT_DATA_CSV: {},{}", var.id(), valuesStr);
            }
        }

//        var.setPonCount(Integer.parseInt(svValues.get(FIELD_INDEX_PON_COUNT)));
//        var.setPonRegionCount(Integer.parseInt(svValues.get(FIELD_INDEX_PON_REGION_COUNT)));
//        var.setLineElements(svValues.get(FIELD_INDEX_LE_START), svValues.get(FIELD_INDEX_LE_END));
//        var.setFragileSites(svValues.get(FIELD_INDEX_FS_START), svValues.get(FIELD_INDEX_FS_END));
    }

}
