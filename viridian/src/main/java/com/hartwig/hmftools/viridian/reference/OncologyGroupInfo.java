package com.hartwig.hmftools.viridian.reference;

import java.util.LinkedHashMap;
import java.util.Map;

import com.hartwig.hmftools.common.utils.file.DelimFileReader;
import com.hartwig.hmftools.common.virus.VirusLikelihoodType;
import com.hartwig.hmftools.common.virus.VirusType;
import com.hartwig.hmftools.viridian.common.UserInputError;

import org.jetbrains.annotations.Nullable;

// Per-oncology-group reference data, mainly the reporting policy data.
public record OncologyGroupInfo(
        OncologyGroup group,
        @Nullable VirusType reportingType,
        @Nullable VirusLikelihoodType driverLikelihood)
{
    public OncologyGroupInfo
    {
        if((reportingType == null) != (driverLikelihood == null))
        {
            throw new IllegalArgumentException("reportingType and driverLikelihood must be set together: " + group);
        }
    }

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
                VirusLikelihoodType driverLikelihood = parseDriverLikelihood(row.getStringOrNull(Columns.driver_likelihood));
                OncologyGroupInfo oncologyGroupInfo = new OncologyGroupInfo(group, reportingType, driverLikelihood);
                OncologyGroupInfo previous = result.put(group, oncologyGroupInfo);
                if(previous != null)
                {
                    throw new UserInputError("Oncology group info has duplicate group: " + group);
                }
            }
        }
        return result;
    }

    @Nullable
    private static VirusLikelihoodType parseDriverLikelihood(@Nullable String value)
    {
        if(value == null)
        {
            return null;
        }
        VirusLikelihoodType likelihood = VirusLikelihoodType.valueOf(value);
        if(likelihood == VirusLikelihoodType.UNKNOWN)
        {
            throw new UserInputError("Oncology group info has invalid driver likelihood: " + value);
        }
        return likelihood;
    }

    private enum Columns
    {
        oncology_group,
        reporting_type,
        driver_likelihood
    }
}
