package com.hartwig.hmftools.svanalysis.annotators;

import static com.hartwig.hmftools.svanalysis.types.SvClusterData.SVI_END;
import static com.hartwig.hmftools.svanalysis.types.SvClusterData.SVI_START;
import static com.hartwig.hmftools.svanalysis.types.SvClusterData.isStart;
import static com.hartwig.hmftools.svanalysis.types.SvGeneData.DRIVER_DEL;
import static com.hartwig.hmftools.svanalysis.types.SvGeneData.DRIVER_TYPE_TSG;

import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.svanalysis.types.SvClusterData;
import com.hartwig.hmftools.svanalysis.types.SvGeneData;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;

public class GeneAnnotator {

    private Map<String, List<SvGeneData>> mSampleGeneData;

//    public static final String KNOWN_FS = "true";

    private static final Logger LOGGER = LogManager.getLogger(GeneAnnotator.class);

    public GeneAnnotator ()
    {
        mSampleGeneData = new HashMap();
    }

    public void loadGeneDriverFile(final String filename)
    {
        if(filename.isEmpty())
            return;

        try {

            BufferedReader fileReader = new BufferedReader(new FileReader(filename));

            String line;

            line = fileReader.readLine(); // skip headers

            int geneDataCount = 0;
            while ((line = fileReader.readLine()) != null) {

                // parse CSV data
                String[] items = line.split(",");

                if(items.length < 20)
                    continue;

                final String sampleId = items[0];

                if(!mSampleGeneData.containsKey(sampleId))
                {
                    mSampleGeneData.put(sampleId, Lists.newArrayList());
                }

                List<SvGeneData> geneDataMap = mSampleGeneData.get(sampleId);

                geneDataMap.add(new SvGeneData(items));
                ++geneDataCount;
            }

            LOGGER.debug("loaded samples({}) and geneData({})", mSampleGeneData.size(), geneDataCount);
        }
        catch(IOException exception)
        {
            LOGGER.error("Failed to read gene data CSV file({})", filename);
        }
    }

    public boolean hasData() { return !mSampleGeneData.isEmpty(); }

    public void addGeneData(final String sampleId, final SvClusterData var)
    {
        if(mSampleGeneData.isEmpty())
            return;

        if(!mSampleGeneData.containsKey(sampleId))
            return;

        final List<SvGeneData> geneDataMap = mSampleGeneData.get(sampleId);

        for(SvGeneData geneData : geneDataMap)
        {
            // for now, only hanle TSG DELs since they ought to match exactly
            if(!geneData.driverType().equals(DRIVER_TYPE_TSG))
                continue;

            if(!geneData.driver().equals(DRIVER_DEL))
                continue;

            // check each breakend if the orientations face away from the deleted region
            for(int be = SVI_START; be <= SVI_END; ++be)
            {
                boolean useStart = isStart(be);

                if(!geneData.chromosome().equals(var.chromosome(useStart)))
                    continue;

//                if((useStart && var.orientation(useStart) != 1) || (!useStart && var.orientation(useStart) != -1))
//                    continue;

                if(geneData.startCNRegion() == var.position(useStart) && var.orientation(useStart) == 1)
                {
                    geneData.addSvData(var, true);
                    var.setGeneData(geneData, useStart);
                }
                else if(geneData.endCNRegion() == var.position(useStart) && var.orientation(useStart) == -1)
                {
                    geneData.addSvData(var, false);
                    var.setGeneData(geneData, useStart);
                }
            }
        }
    }

    public void reportGeneMatchData(final String sampleId)
    {
        if(mSampleGeneData.isEmpty())
            return;

        if(!mSampleGeneData.containsKey(sampleId))
            return;

        final List<SvGeneData> geneDataMap = mSampleGeneData.get(sampleId);

        for(SvGeneData geneData : geneDataMap) {

            if (!geneData.driverType().equals(DRIVER_TYPE_TSG))
                continue;

            if(!geneData.driver().equals(DRIVER_DEL))
                continue;

            if(geneData.getStartSvList().isEmpty()) {
                LOGGER.debug("sample({}) gene({}) start unmatched with regionSupport({})", sampleId, geneData.gene(), geneData.startRegionType());
            }

            if(geneData.getEndSvList().isEmpty())
            {
                LOGGER.debug("sample({}) gene({}) end unmatched with regionSupport({})", sampleId, geneData.gene(), geneData.endRegionType());
            }

            if(!geneData.getStartSvList().isEmpty() && !geneData.getEndSvList().isEmpty())
            {
                if(geneData.getStartSvList().size() == geneData.getEndSvList().size())
                {
                    if(geneData.getStartSvList().get(0).equals(geneData.getEndSvList().get(0))) {
                        LOGGER.info("sample({}) gene({}) matches single SV({})", sampleId, geneData.gene(), geneData.getEndSvList().get(0).id());
                    }
                    else {
                        LOGGER.info("sample({}) gene({}) matches diff SVs({} & {})",
                                sampleId, geneData.gene(), geneData.getStartSvList().get(0).id(), geneData.getEndSvList().get(0).id());
                    }
                }
                else
                {
                    LOGGER.info("sample({}) gene({}) matches multiple SVs, counts({} & {})",
                            geneData.getStartSvList().get(0).id(), geneData.getStartSvList().size(), geneData.getEndSvList().size());
                }
            }
        }
    }

}
