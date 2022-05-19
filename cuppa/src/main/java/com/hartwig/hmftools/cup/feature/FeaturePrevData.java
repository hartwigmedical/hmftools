package com.hartwig.hmftools.cup.feature;

import static java.lang.String.format;

import static com.hartwig.hmftools.cup.CuppaConfig.DATA_DELIM;

public class FeaturePrevData
{
    public final String CancerType;
    public final String Name;
    public final FeatureType Type;
    public final double RawPrevalence;
    public double Prevalence;

    public FeaturePrevData(final String cancerType, final String name, final FeatureType type, final double prevalence)
    {
        CancerType = cancerType;
        Name = name;
        Type = type;
        RawPrevalence = prevalence;
        Prevalence = prevalence;
    }

    public static FeaturePrevData from(final String data)
    {
        final String[] items = data.split(DATA_DELIM, -1);
        if(items.length != 4)
            return null;

        return new FeaturePrevData(items[0], items[1], FeatureType.valueOf(items[2]), Double.parseDouble(items[3]));
    }

    public String toString() { return format("ct(%s) gene(%s) type(%s) prev(%.4f adj=%.4f)",
            CancerType, Name, Type, RawPrevalence, Prevalence); }

    public String typeName() { return featureTypeName(Type, Name); }

    public static String TYPE_NAME_DELIM = "_";

    public static String featureTypeName(final FeatureType type, final String name) { return type + TYPE_NAME_DELIM + name; }

    public static FeatureType typeFromFeatureTypeName(final String typeName)
    {
        return FeatureType.valueOf(typeName.substring(0, typeName.indexOf(TYPE_NAME_DELIM)));
    }

    public static String nameFromfeatureTypeName(final String typeName)
    {
        return typeName.substring(typeName.indexOf(TYPE_NAME_DELIM) + 1);
    }
}
