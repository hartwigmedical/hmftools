package com.hartwig.hmftools.peach.effect;

import java.util.List;
import java.util.Set;
import java.util.stream.Collectors;

import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

public class DrugInfoStore
{
    @NotNull
    private final List<DrugInfo> drugInfos;

    public DrugInfoStore(@NotNull List<DrugInfo> drugInfos)
    {
        this.drugInfos = drugInfos;
    }

    @NotNull
    public Set<String> getRelevantDrugNames(@NotNull String geneName)
    {
        return drugInfos.stream().filter(i -> i.geneName().equals(geneName)).map(DrugInfo::drugName).collect(Collectors.toSet());
    }

    @Nullable
    public String getPrescriptionInfoUrl(@NotNull String geneName, @NotNull String drugName)
    {
        DrugInfo matchingDrugInfo = getMatchingDrugInfo(geneName, drugName);
        return matchingDrugInfo == null ? null : matchingDrugInfo.prescriptionInfoUrl();
    }

    @Nullable
    private DrugInfo getMatchingDrugInfo(@NotNull String geneName, @NotNull String drugName)
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
