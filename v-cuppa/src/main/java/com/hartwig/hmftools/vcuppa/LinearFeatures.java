package com.hartwig.hmftools.vcuppa;

import java.util.LinkedHashMap;
import java.util.Map;

// linear input features. They are drivers and signatures.
public class LinearFeatures
{
    // this value must match the value used in vcuppa_train.py
    public static double DRIVER_FILL_NAN = -0.2;

    public static class Feature
    {
        private final String key;
        private final boolean log1p;
        private double value = DRIVER_FILL_NAN;

        public Feature(String key, boolean log1p)
        {
            this.key = key;
            this.log1p = log1p;
        }

        public String getKey()
        {
            return key;
        }

        public boolean isLog1p()
        {
            return log1p;
        }

        public double getValue()
        {
            return value;
        }

        public void setValue(double value)
        {
            this.value = value;
        }
    }

    // we want the values to be in insertion order
    private final Map<String, Feature> features = new LinkedHashMap<>();

    public LinearFeatures()
    {
    }

    public void addFeatureKey(String key, boolean log1p)
    {
        features.put(key, new Feature(key, log1p));
    }

    public void setFeatureValueIfPresent(String key, double value)
    {
        // apply log1p here
        Feature f = features.get(key);
        if (f != null)
        {
            f.setValue(f.log1p ? Math.log1p(value) : value);
        }
    }

    public float[] toFloatArray()
    {
        float[] array = new float[features.size()];
        int index = 0;
        for(Feature feature : features.values())
        {
            array[index++] = (float)feature.getValue();
        }
        return array;
    }
}
