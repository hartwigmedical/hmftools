package com.hartwig.hmftools.viridian.reference;

import java.util.LinkedHashMap;
import java.util.Map;

import com.hartwig.hmftools.common.utils.file.DelimFileReader;
import com.hartwig.hmftools.common.virus.VirusType;
import com.hartwig.hmftools.viridian.common.UserInputError;

import org.jetbrains.annotations.Nullable;

// Per-virus-group reference data, mainly the reporting policy data.
public record VirusGroupInfo(
        VirusGroup group,
        @Nullable VirusType reportingType
)
{
    public boolean isReportable()
    {
        return reportingType != null;
    }

    public static Map<VirusGroup, VirusGroupInfo> load(String tsvFile)
    {
        Map<VirusGroup, VirusGroupInfo> result = new LinkedHashMap<>();
        try(DelimFileReader reader = new DelimFileReader(tsvFile))
        {
            for(DelimFileReader.Row row : reader)
            {
                VirusGroup group = new VirusGroup(row.getString(Columns.virus_group));
                String reportingTypeValue = row.getStringOrNull(Columns.reporting_type);
                VirusType reportingType = reportingTypeValue == null ? null : VirusType.fromVirusName(reportingTypeValue);
                VirusGroupInfo virusGroupInfo = new VirusGroupInfo(group, reportingType);
                VirusGroupInfo previous = result.put(group, virusGroupInfo);
                if(previous != null)
                {
                    throw new UserInputError("Virus group info has duplicate group: " + group);
                }
            }
        }
        return result;
    }

    private enum Columns
    {
        virus_group,
        reporting_type
    }
}
