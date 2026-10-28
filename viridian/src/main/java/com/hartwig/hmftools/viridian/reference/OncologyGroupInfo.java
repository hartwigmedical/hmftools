package com.hartwig.hmftools.viridian.reference;

import java.util.LinkedHashMap;
import java.util.Map;

import com.hartwig.hmftools.common.utils.file.DelimFileReader;
import com.hartwig.hmftools.common.virus.VirusType;
import com.hartwig.hmftools.viridian.common.UserInputError;

import org.jetbrains.annotations.Nullable;

// Per-oncology-group reference data, mainly the reporting policy data.
public record OncologyGroupInfo(
        OncologyGroup group,
        @Nullable VirusType reportingType
)
{
    public boolean isReportable()
    {
        return reportingType != null;
    }

    public static Map<OncologyGroup, OncologyGroupInfo> load(String tsvFile)
    {
        Map<OncologyGroup, OncologyGroupInfo> result = new LinkedHashMap<>();
        try(DelimFileReader reader = new DelimFileReader(tsvFile))
        {
            for(DelimFileReader.Row row : reader)
            {
                OncologyGroup group = new OncologyGroup(row.getString(Columns.oncology_group));
                String reportingTypeValue = row.getStringOrNull(Columns.reporting_type);
                VirusType reportingType = reportingTypeValue == null ? null : VirusType.fromVirusName(reportingTypeValue);
                OncologyGroupInfo oncologyGroupInfo = new OncologyGroupInfo(group, reportingType);
                OncologyGroupInfo previous = result.put(group, oncologyGroupInfo);
                if(previous != null)
                {
                    throw new UserInputError("Oncology group info has duplicate group: " + group);
                }
            }
        }
        return result;
    }

    private enum Columns
    {
        oncology_group,
        reporting_type
    }
}
