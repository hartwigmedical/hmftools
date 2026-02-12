package com.hartwig.hmftools.datamodel.finding;

public class FindingsJson extends JsonReadWriter<FindingRecord>
{
    public FindingsJson()
    {
        super(FindingRecord.class);
    }
}
