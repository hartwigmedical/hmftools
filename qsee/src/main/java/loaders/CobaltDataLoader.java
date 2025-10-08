package loaders;

import static common.QSeeConstants.QC_LOGGER;

import java.io.IOException;
import java.util.ArrayList;
import java.util.List;
import java.util.Map;

import com.hartwig.hmftools.common.cobalt.CobaltGcMedianFile;
import com.hartwig.hmftools.common.cobalt.GcMedianReadDepth;
import com.hartwig.hmftools.common.genome.gc.GCBucket;

import prep.DataType;
import prep.FeatureType;
import prep.FeatureValue;
import prep.PrepConfig;

public class CobaltDataLoader
{
    PrepConfig mConfig;

    public CobaltDataLoader(PrepConfig config)
    {
        mConfig = config;
    }

    public List<FeatureValue> loadCobaltGcMedians(String sampleId)
    {

        try {
            String filePath = CobaltGcMedianFile.generateFilename(mConfig.getCobaltDir(sampleId), sampleId);

            GcMedianReadDepth gcMedianReadDepth = CobaltGcMedianFile.read(filePath);

            Map<GCBucket, Double> medianReadDepths = gcMedianReadDepth.medianReadDepthPerGCBucket();
            double overallMedianReadDepth = gcMedianReadDepth.medianReadDepth();

            List<FeatureValue> featureValues = new ArrayList<>();

            for(GCBucket bucket : medianReadDepths.keySet())
            {
                double medianReadDepth = medianReadDepths.get(bucket);

                if(medianReadDepth == GcMedianReadDepth.NO_READ_DEPTH_VALUE)
                    medianReadDepth = 0;

                double normalisedDepth = medianReadDepth / overallMedianReadDepth;

                FeatureValue featureValue = new FeatureValue(
                        bucket.toString(),
                        Double.toString(normalisedDepth),
                        DataType.NUMBER,
                        FeatureType.COBALT_GC_MEDIAN
                );

                featureValues.add(featureValue);
            }

            return featureValues;
        }
        catch (IOException e)
        {
            QC_LOGGER.error("Failed to load Cobalt GC median read depth file: {}", e.toString());
            e.printStackTrace();
            return null;
        }
    }
}
