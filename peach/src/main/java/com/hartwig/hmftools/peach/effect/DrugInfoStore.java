package com.hartwig.hmftools.peach.effect;

import java.util.List;
import java.util.Set;
import java.util.stream.Collectors;

import org.jetbrains.annotations.NotNull;

public class DrugInfoStore
{
    @NotNull
    private final List<DrugInfo> drugInfos;

    public DrugInfoStore(@NotNull List<DrugInfo> drugInfos)
    {
        this.drugInfos = drugInfos;
    }

    public Set<String> getRelevantDrugNames(@NotNull String geneName)
    {
        return drugInfos.stream().filter(i -> i.geneName().equals(geneName)).map(DrugInfo::drugName).collect(Collectors.toSet());
    }

    public Set<String> getGeneralInfoUrls(@NotNull String geneName, @NotNull String drugName)
    {
        return drugInfos.stream()
                .filter(i -> i.geneName().equals(geneName) && i.drugName().equals(drugName))
                .map(DrugInfo::generalInfoUrl)
                .collect(Collectors.toSet());
    }

    public Set<String> getPrescriptionInfoUrls(@NotNull String geneName, @NotNull String drugName)
    {
        return drugInfos.stream()
                .filter(i -> i.geneName().equals(geneName) && i.drugName().equals(drugName))
                .map(DrugInfo::prescriptionInfoUrl)
                .collect(Collectors.toSet());
    }
}
