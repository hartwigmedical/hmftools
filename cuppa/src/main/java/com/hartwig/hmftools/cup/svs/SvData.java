package com.hartwig.hmftools.cup.svs;

import static com.hartwig.hmftools.cup.SampleAnalyserConfig.DATA_DELIM;
import static com.hartwig.hmftools.cup.svs.SvDataType.DUP_LONG;
import static com.hartwig.hmftools.cup.svs.SvDataType.DUP_MEDIUM;
import static com.hartwig.hmftools.cup.svs.SvDataType.DUP_SHORT;
import static com.hartwig.hmftools.cup.svs.SvDataType.FRAGILE_SITE;
import static com.hartwig.hmftools.cup.svs.SvDataType.LINE;
import static com.hartwig.hmftools.cup.svs.SvDataType.MAX_EVENT_SIZE;
import static com.hartwig.hmftools.cup.svs.SvDataType.TELOMERIC_SGL;
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

        // SampleId,Line,FragileSite,DupShort,DupMedium,DupLong,MaxEventSize,TelomericSgl

        SvData svData = new SvData(items[fieldsIndexMap.get("SampleId")]);
        svData.setCount(LINE, Integer.parseInt(items[fieldsIndexMap.get("Line")]));
        svData.setCount(FRAGILE_SITE, Integer.parseInt(items[fieldsIndexMap.get("FragileSite")]));
        svData.setCount(DUP_SHORT, Integer.parseInt(items[fieldsIndexMap.get("DupShort")]));
        svData.setCount(DUP_MEDIUM, Integer.parseInt(items[fieldsIndexMap.get("DupMedium")]));
        svData.setCount(DUP_LONG, Integer.parseInt(items[fieldsIndexMap.get("DupLong")]));
        svData.setCount(MAX_EVENT_SIZE, Integer.parseInt(items[fieldsIndexMap.get("MaxEventSize")]));
        svData.setCount(TELOMERIC_SGL, Integer.parseInt(items[fieldsIndexMap.get("TelomericSgl")]));

        return svData;

    }

}
