package com.hartwig.hmftools.cup.feature;

import static com.hartwig.hmftools.cup.CuppaConfig.DATA_DELIM;
import static com.hartwig.hmftools.cup.CuppaConfig.FLD_SAMPLE_ID;
import static com.hartwig.hmftools.cup.CuppaConfig.SUBSET_DELIM;
import static com.hartwig.hmftools.cup.feature.FeaturePrevData.featureTypeName;

import java.util.Map;
import java.util.StringJoiner;

import com.google.common.collect.Maps;

public class SampleFeatureData
{
    public final String SampleId;
    public final String Name; // driver gene, fusion genes or virus
    public final FeatureType Type;
    public final double Likelihood;
    public final Map<String,String> ExtraInfo;

    public static final String DRIVER_TYPE = "DriverType";
    public static final String DRIVER_TYPE_AMP = "AMP";
    public static final String DRIVER_TYPE_DEL = "DEL";

    public static final String DRIVER_CHROMOSOME = "DriverChromosome";

    public SampleFeatureData(final String sampleId, final String name, final FeatureType type, final double likelihood)
    {
        SampleId = sampleId;
        Name = name;
        Type = type;
        Likelihood = likelihood;
        ExtraInfo = Maps.newHashMap();
    }

    public String typeName() { return featureTypeName(Type, Name); }

    public static String header()
    {
        StringJoiner sj = new StringJoiner(DATA_DELIM);
        sj.add(FLD_SAMPLE_ID);
        sj.add("Name");
        sj.add("Type");
        sj.add("Likelihood");
        sj.add("ExtraInfo");
        return sj.toString();
    }

    public String toCsv()
    {
        StringJoiner sj = new StringJoiner(DATA_DELIM);
        sj.add(Name);
        sj.add(Type.toString());
        sj.add(String.valueOf(Likelihood));

        StringJoiner sjEI = new StringJoiner(SUBSET_DELIM);

        for(Map.Entry<String,String> entry : ExtraInfo.entrySet())
        {
            sjEI.add(String.format("%s=%s", entry.getKey(), entry.getValue()));
        }

        sj.add(sjEI.toString());
        return sj.toString();
    }

    public static SampleFeatureData from(final String data)
    {
        final String[] items = data.split(DATA_DELIM, -1);
        if(items.length != 5)
            return null;

        SampleFeatureData driverData = new SampleFeatureData(items[0], items[1], FeatureType.valueOf(items[2]), Double.parseDouble(items[3]));

        final String extraInfoStr = items[4];
        if(!extraInfoStr.isEmpty())
        {
            String[] eiItems = extraInfoStr.split(SUBSET_DELIM, -1);
            for(String eiItem : eiItems)
            {
                String[] eiInfo = eiItem.split("=");
                driverData.ExtraInfo.put(eiInfo[0], eiInfo[1]);
            }
        }

        return driverData;
    }

    public String toString() { return String.format("sample(%s) gene(%s) type(%s)", SampleId, Name, Type); }

    public String getExtraInfo(final String key) { return ExtraInfo.containsKey(key) ? ExtraInfo.get(key) : ""; }
}
