package com.hartwig.hmftools.cup.svs;

import static com.hartwig.hmftools.cup.CuppaConfig.DATA_DELIM;
import static com.hartwig.hmftools.cup.CuppaConfig.SAMPLE_ID;
import static com.hartwig.hmftools.cup.svs.SvDataType.typeIndex;

import java.util.Arrays;
import java.util.Map;
import java.util.StringJoiner;

public class SvData
{
    public final String SampleId;
    public final int[] TypeCounts;

    public SvData(final String sampleId)
    {
        SampleId = sampleId;
        TypeCounts = new int[SvDataType.count()];
    }

    public int getCount(final SvDataType type) { return TypeCounts[typeIndex(type)]; }
    public void setCount(final SvDataType type, int count) { TypeCounts[typeIndex(type)] = count; }

    public static SvData from(final Map<String,Integer> fieldsIndexMap, final String data)
    {
        final String[] items = data.split(DATA_DELIM, -1);

        SvData svData = new SvData(items[fieldsIndexMap.get(SAMPLE_ID)]);

        for(SvDataType type : SvDataType.values())
        {
            svData.setCount(type, Integer.parseInt(items[fieldsIndexMap.get(type.toString())]));
        }

        return svData;
    }

    public static String header()
    {
        StringJoiner sj = new StringJoiner(DATA_DELIM);
        sj.add(SAMPLE_ID);

        for(SvDataType type : SvDataType.values())
        {
            sj.add(type.toString());
        }

        return sj.toString();
    }

    public String toCsv()
    {
        StringJoiner sj = new StringJoiner(DATA_DELIM);

        for(int i = 0; i < TypeCounts.length; ++i)
        {
            sj.add(String.valueOf(TypeCounts[i]));
        }

        return sj.toString();
    }

}
