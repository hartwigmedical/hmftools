package com.hartwig.hmftools.finding;

import java.io.IOException;
import java.nio.file.Path;
import java.util.Map;
import java.util.function.Function;
import java.util.stream.Collectors;

import com.hartwig.hmftools.common.driver.panel.DriverGene;
import com.hartwig.hmftools.common.driver.panel.DriverGeneFile;
import com.hartwig.hmftools.common.purple.Gender;
import com.hartwig.hmftools.datamodel.orange.OrangeRefGenomeVersion;
import com.hartwig.hmftools.finding.clinicaltranscript.ClinicalTranscriptFile;
import com.hartwig.hmftools.finding.clinicaltranscript.ClinicalTranscriptsModel;

import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

public class FindingConfig
{
    @Nullable
    private final ClinicalTranscriptsModel clinicalTranscriptsModel;
    @NotNull
    private final Map<String, DriverGene> driverGenes;
    @Nullable
    private final Gender gender;

    public static FindingConfig createFindingConfig(@Nullable Path clinicalTranscriptsTsv,
            @Nullable Path driverGeneTsv, OrangeRefGenomeVersion orangeRefGenomeVersion,
            @Nullable Gender gender) throws IOException
    {
        ClinicalTranscriptsModel clinicalTranscriptsModel = clinicalTranscriptsTsv != null ?
                ClinicalTranscriptFile.buildFromTsv(orangeRefGenomeVersion, clinicalTranscriptsTsv) : null;
        Map<String, DriverGene> driverGenes = driverGenesMap(driverGeneTsv);
        return new FindingConfig(clinicalTranscriptsModel, driverGenes, gender);
    }

    private static Map<String, DriverGene> driverGenesMap(@Nullable Path driverGeneTsv) throws IOException
    {
        return driverGeneTsv != null ? DriverGeneFile.read(driverGeneTsv)
                .stream()
                .collect(Collectors.toMap(DriverGene::gene, Function.identity())) : Map.of();
    }

    public FindingConfig(@Nullable final ClinicalTranscriptsModel clinicalTranscriptsModel, final Map<String, DriverGene> driverGenes,
            @Nullable Gender gender)
    {
        this.clinicalTranscriptsModel = clinicalTranscriptsModel;
        this.driverGenes = driverGenes;
        this.gender = gender;
    }

    @Nullable
    public String findCanonicalTranscriptForGene(String gene)
    {
        return clinicalTranscriptsModel != null ? clinicalTranscriptsModel.findCanonicalTranscriptForGene(gene) : null;
    }

    @Nullable
    public DriverGene getDriverGene(String gene)
    {
        return driverGenes.get(gene);
    }

    public Gender gender()
    {
        return gender;
    }
}
