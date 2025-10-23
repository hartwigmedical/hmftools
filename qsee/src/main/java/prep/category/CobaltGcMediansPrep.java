package prep.category;

import static common.QSeeConstants.QC_LOGGER;

import java.io.IOException;
import java.util.ArrayList;
import java.util.Comparator;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;

import com.hartwig.hmftools.common.cobalt.CobaltGcMedianFile;
import com.hartwig.hmftools.common.cobalt.GcMedianReadDepth;
import com.hartwig.hmftools.common.genome.gc.GCBucket;

import org.apache.commons.lang3.tuple.Pair;

import feature.FeatureKey;
import feature.FeatureType;
import feature.Feature;
import prep.CategoryPrep;
import prep.PrepConfig;

public class CobaltGcMediansPrep implements CategoryPrep
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

    public static List<Feature> normaliseMedianReadDepths(GcMedianReadDepth gcMedianReadDepth)
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

            FeatureKey key = FeatureKey.of(
                    FeatureType.COBALT_GC_MEDIAN,
                    Pair.of(KEY_FLD_GC_BUCKET, String.valueOf(bucket.bucket()))
            );

            Feature feature = new Feature(key, normalisedDepth);

            features.add(feature);
        }

        return features;
    }

    public List<Feature> extractSampleData(String sampleId)
    {
        GcMedianReadDepth gcMedianReadDepth = loadCobaltGcMedianFile(sampleId);
        List<Feature> features = normaliseMedianReadDepths(gcMedianReadDepth);

        return features;
    }
}
