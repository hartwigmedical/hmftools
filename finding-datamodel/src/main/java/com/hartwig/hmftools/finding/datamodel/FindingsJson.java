package com.hartwig.hmftools.finding.datamodel;

public class FindingsJson extends JsonReadWriter<FindingRecord>
{
    public FindingsJson()
    {
        super(FindingRecord.class);
    }
}
