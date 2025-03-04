package com.hartwig.hmftools.peach.effect;

import java.util.List;
import java.util.Set;
import java.util.stream.Collectors;

import org.jetbrains.annotations.Nullable;

public class DrugInfoStore
{
    private final List<DrugInfo> drugInfos;

    public DrugInfoStore(final List<DrugInfo> drugInfos)
    {
        this.drugInfos = drugInfos;
    }

    public Set<String> getRelevantDrugNames(final String geneName)
    {
        return drugInfos.stream().filter(i -> i.geneName().equals(geneName)).map(DrugInfo::drugName).collect(Collectors.toSet());
    }

    @Nullable
    public String getPrescriptionInfoUrl(final String geneName, final String drugName)
    {
        DrugInfo matchingDrugInfo = getMatchingDrugInfo(geneName, drugName);
        return matchingDrugInfo == null ? null : matchingDrugInfo.prescriptionInfoUrl();
    }

    @Nullable
    private DrugInfo getMatchingDrugInfo(final String geneName, final String drugName)
    {
        List<DrugInfo> matchingDrugInfos =
                drugInfos.stream().filter(i -> i.geneName().equals(geneName) && i.drugName().equals(drugName)).collect(Collectors.toList());
        if(matchingDrugInfos.isEmpty())
        {
            return null;
        }
        else if(matchingDrugInfos.size() == 1)
        {
            return matchingDrugInfos.get(0);
        }
        else
        {
            String errorMessage = String.format("Multiple urls configured for drug %s with gene %s", drugName, geneName);
            throw new IllegalStateException(errorMessage);
        }
    }
}
