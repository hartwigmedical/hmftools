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
import com.hartwig.hmftools.finding.clinicalrelevantgenecopynumber.ClinicalRelevantGeneCopyNumberFile;
import com.hartwig.hmftools.finding.clinicalrelevantgenecopynumber.ClinicalRelevantGeneCopyNumberModel;
import com.hartwig.hmftools.finding.clinicaltranscript.ClinicalTranscriptFile;
import com.hartwig.hmftools.finding.clinicaltranscript.ClinicalTranscriptsModel;

import org.jetbrains.annotations.Nullable;

record FindingConfig(@Nullable ClinicalTranscriptsModel clinicalTranscriptsModel,
                     @Nullable ClinicalRelevantGeneCopyNumberModel clinicalRelevantGeneCopyNumberModel,
                     Map<String, DriverGene> driverGenes,
                     @Nullable Gender gender,
                     boolean geneCopyNumbersOptional)
{
    public static FindingConfig createFindingConfig(@Nullable Path clinicalTranscriptsTsv,
            @Nullable Path clinicalRelevantGeneCopyNumbersTsv,
            @Nullable Path driverGeneTsv, OrangeRefGenomeVersion orangeRefGenomeVersion,
            @Nullable Gender gender, boolean geneCopyNumbersOptional) throws IOException
    {
        ClinicalTranscriptsModel clinicalTranscriptsModel = clinicalTranscriptsTsv != null ?
                ClinicalTranscriptFile.buildFromTsv(orangeRefGenomeVersion, clinicalTranscriptsTsv) : null;
        ClinicalRelevantGeneCopyNumberModel clinicalRelevantGeneCopyNumberModel = clinicalRelevantGeneCopyNumbersTsv != null ?
                ClinicalRelevantGeneCopyNumberFile.buildFromTsv(clinicalRelevantGeneCopyNumbersTsv) : null;
        Map<String, DriverGene> driverGenes = driverGenesMap(driverGeneTsv);
        return new FindingConfig(clinicalTranscriptsModel, clinicalRelevantGeneCopyNumberModel, driverGenes, gender, geneCopyNumbersOptional);
    }

    private static Map<String, DriverGene> driverGenesMap(@Nullable Path driverGeneTsv) throws IOException
    {
        return driverGeneTsv != null ? DriverGeneFile.read(driverGeneTsv)
                .stream()
                .collect(Collectors.toMap(DriverGene::gene, Function.identity())) : Map.of();
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

    public boolean findClinicalRelevantGeneCopyNumber(String gene)
    {
        return clinicalRelevantGeneCopyNumberModel != null && clinicalRelevantGeneCopyNumberModel.findClinicalRelevantCopyNumberGene(gene);
    }
}