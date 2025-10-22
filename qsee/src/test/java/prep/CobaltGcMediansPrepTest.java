package prep;

import static org.junit.Assert.assertEquals;

import java.util.List;
import java.util.Map;

import com.hartwig.hmftools.common.cobalt.GcMedianReadDepth;
import com.hartwig.hmftools.common.genome.gc.ImmutableGCBucket;

import org.junit.Test;

import feature.Feature;
import prep.category.CobaltGcMediansPrep;

public class CobaltGcMediansPrepTest
{
    @Test
    public void canNormaliseMedianReadDepths(){
        GcMedianReadDepth gcMedianReadDepth = new GcMedianReadDepth(
                100.0, 100.0,
                Map.of(new ImmutableGCBucket(25), 100.0,
                       new ImmutableGCBucket(26), 105.0)
        );

        List<Feature> features = CobaltGcMediansPrep.normaliseMedianReadDepths(gcMedianReadDepth);

        Feature actualFeature;

        actualFeature = features.get(0);
        assertEquals("GCBucket=25", actualFeature.mKey);
        assertEquals(1.0, actualFeature.mValue, 0.001);

        actualFeature = features.get(1);
        assertEquals("GCBucket=26", actualFeature.mKey);
        assertEquals(1.05, actualFeature.mValue, 0.001);
    }
}
