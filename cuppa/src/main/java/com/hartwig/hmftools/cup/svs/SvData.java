package com.hartwig.hmftools.cup.svs;

import static com.hartwig.hmftools.cup.SampleAnalyserConfig.DATA_DELIM;
import static com.hartwig.hmftools.cup.svs.SvDataType.typeIndex;

import java.util.Map;

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

        SvData svData = new SvData(items[fieldsIndexMap.get("SampleId")]);

        for(SvDataType type : SvDataType.values())
        {
            svData.setCount(type, Integer.parseInt(items[fieldsIndexMap.get(type.toString())]));
        }

        return svData;

    }

}
