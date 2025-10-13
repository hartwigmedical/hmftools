package loaders;

import static common.QSeeConstants.QC_LOGGER;

import java.io.IOException;
import java.util.ArrayList;
import java.util.List;
import java.util.Map;

import com.hartwig.hmftools.common.cobalt.CobaltGcMedianFile;
import com.hartwig.hmftools.common.cobalt.GcMedianReadDepth;
import com.hartwig.hmftools.common.genome.gc.GCBucket;

import prep.FeatureType;
import prep.FeatureValue;
import prep.PrepConfig;

public class CobaltGcMediansPrep
{
    PrepConfig mConfig;

    public CobaltGcMediansPrep(PrepConfig config)
    {
        mConfig = config;
    }

    public List<FeatureValue<Double>> extractSampleData(String sampleId)
    {
        try {
            String filePath = CobaltGcMedianFile.generateFilename(mConfig.getCobaltDir(sampleId), sampleId);

            GcMedianReadDepth gcMedianReadDepth = CobaltGcMedianFile.read(filePath);

            Map<GCBucket, Double> medianReadDepths = gcMedianReadDepth.medianReadDepthPerGCBucket();
            double overallMedianReadDepth = gcMedianReadDepth.medianReadDepth();

            List<FeatureValue<Double>> featureValues = new ArrayList<>();

            for(GCBucket bucket : medianReadDepths.keySet())
            {
                double medianReadDepth = medianReadDepths.get(bucket);

                if(medianReadDepth == GcMedianReadDepth.NO_READ_DEPTH_VALUE)
                    medianReadDepth = 0;

                double normalisedDepth = medianReadDepth / overallMedianReadDepth;

                FeatureValue<Double> featureValue = new FeatureValue<>(
                        bucket.toString(),
                        normalisedDepth,
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
