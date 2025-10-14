package prep;

import static common.QSeeConstants.QC_LOGGER;

import java.io.IOException;
import java.util.ArrayList;
import java.util.List;
import java.util.Map;

import com.hartwig.hmftools.common.cobalt.CobaltGcMedianFile;
import com.hartwig.hmftools.common.cobalt.GcMedianReadDepth;
import com.hartwig.hmftools.common.genome.gc.GCBucket;

import feature.FeatureType;
import feature.FeatureValue;

public class CobaltGcMediansPrep implements CategoryPrep<Double>
{
    PrepConfig mConfig;

    private static final String KEY_FLD_GC_BUCKET = "GCBucket";

    public CobaltGcMediansPrep(PrepConfig config)
    {
        mConfig = config;
    }

    private GcMedianReadDepth loadCobaltGcMedianFile(String sampleId)
    {
        try
        {
            String filePath = CobaltGcMedianFile.generateFilename(mConfig.getCobaltDir(sampleId), sampleId);
            return CobaltGcMedianFile.read(filePath);
        }
        catch(IOException e)
        {
            QC_LOGGER.error("Failed to load Cobalt GC median read depth file: {}", e.toString());
            e.printStackTrace();
            return null;
        }
    }

    public static List<FeatureValue<Double>> normaliseMedianReadDepths(GcMedianReadDepth gcMedianReadDepth)
    {
        Map<GCBucket, Double> medianReadDepths = gcMedianReadDepth.medianReadDepthPerGCBucket();

        double overallMedianReadDepth = gcMedianReadDepth.medianReadDepth();

        List<FeatureValue<Double>> featureValues = new ArrayList<>();

        for(GCBucket bucket : medianReadDepths.keySet())
        {
            double medianReadDepth = medianReadDepths.get(bucket);
            if(medianReadDepth == GcMedianReadDepth.NO_READ_DEPTH_VALUE)
            {
                medianReadDepth = 0;
            }

            double normalisedDepth = medianReadDepth / overallMedianReadDepth;

            FeatureValue<Double> featureValue = new FeatureValue<>(
                    FeatureValue.keyFromPair(KEY_FLD_GC_BUCKET, String.valueOf(bucket.bucket())),
                    normalisedDepth,
                    FeatureType.COBALT_GC_MEDIAN
            );

            featureValues.add(featureValue);
        }

        return featureValues;
    }

    public List<FeatureValue<Double>> extractSampleData(String sampleId)
    {
        GcMedianReadDepth gcMedianReadDepth = loadCobaltGcMedianFile(sampleId);
        List<FeatureValue<Double>> featureValues = normaliseMedianReadDepths(gcMedianReadDepth);

        return featureValues;
    }
}
