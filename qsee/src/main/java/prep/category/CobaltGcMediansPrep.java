package prep.category;

import static common.QseeConstants.QC_LOGGER;

import java.io.IOException;
import java.util.ArrayList;
import java.util.Comparator;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;

import com.google.common.annotations.VisibleForTesting;
import com.hartwig.hmftools.common.cobalt.CobaltGcMedianFile;
import com.hartwig.hmftools.common.cobalt.GcMedianReadDepth;
import com.hartwig.hmftools.common.genome.gc.GCBucket;

import feature.FeatureKey;
import feature.FeatureType;
import feature.Feature;
import feature.SourceTool;
import prep.CategoryPrep;
import prep.PrepConfig;

public class CobaltGcMediansPrep implements CategoryPrep
{
    private final PrepConfig mConfig;

    private static final String FIELD_GC_BUCKET = "GCBucket";

    public CobaltGcMediansPrep(PrepConfig config)
    {
        mConfig = config;
    }

    private GcMedianReadDepth loadCobaltGcMedianFile(String sampleId) throws IOException
    {
        try
        {
            String filePath = CobaltGcMedianFile.generateFilename(mConfig.getCobaltDir(sampleId), sampleId);
            return CobaltGcMedianFile.read(filePath);
        }
        catch(IOException e)
        {
            QC_LOGGER.error("Failed to load Cobalt GC median read depth file for sample: {}", sampleId, e);
            throw e;
        }
    }

    @VisibleForTesting
    static List<Feature> normaliseMedianReadDepths(GcMedianReadDepth gcMedianReadDepth)
    {
        double overallMedianReadDepth = gcMedianReadDepth.medianReadDepth();

        Map<GCBucket, Double> medianReadDepths = gcMedianReadDepth.medianReadDepthPerGCBucket();

        LinkedHashMap<GCBucket, Double> orderedMedianReadDepths = medianReadDepths.entrySet().stream()
                .sorted(Comparator.comparingInt(entry -> entry.getKey().bucket()))
                .collect(LinkedHashMap::new,
                        (map, entry) -> map.put(entry.getKey(), entry.getValue()),
                        LinkedHashMap::putAll);

        List<Feature> features = new ArrayList<>();

        for(GCBucket bucket : orderedMedianReadDepths.keySet())
        {
            double medianReadDepth = orderedMedianReadDepths.get(bucket);
            if(medianReadDepth == GcMedianReadDepth.NO_READ_DEPTH_VALUE)
            {
                medianReadDepth = 0;
            }

            double normalisedDepth = medianReadDepth / overallMedianReadDepth;

            String featureName = FeatureKey.formMultiFieldName(FIELD_GC_BUCKET, String.valueOf(bucket.bucket()));
            Feature feature = new Feature(FeatureType.GC_BIAS, featureName, normalisedDepth, SourceTool.COBALT);

            features.add(feature);
        }

        return features;
    }

    public List<Feature> extractSampleData(String sampleId) throws IOException
    {
        GcMedianReadDepth gcMedianReadDepth = loadCobaltGcMedianFile(sampleId);
        return normaliseMedianReadDepths(gcMedianReadDepth);
    }
}
