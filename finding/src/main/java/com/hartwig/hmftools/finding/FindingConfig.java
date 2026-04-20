package com.hartwig.hmftools.finding;

import java.io.IOException;
import java.nio.file.Path;
import java.util.Map;
import java.util.function.Function;
import java.util.stream.Collectors;

import com.hartwig.hmftools.common.driver.panel.DriverGene;
import com.hartwig.hmftools.common.driver.panel.DriverGeneFile;
import com.hartwig.hmftools.datamodel.orange.OrangeRefGenomeVersion;
import com.hartwig.hmftools.finding.clinicaltranscript.ClinicalTranscriptFile;
import com.hartwig.hmftools.finding.clinicaltranscript.ClinicalTranscriptsModel;

import org.jetbrains.annotations.Nullable;

public class FindingConfig
{
    @Nullable
    private final ClinicalTranscriptsModel clinicalTranscriptsModel;

    public static FindingConfig createFindingConfig(@Nullable Path clinicalTranscriptsTsv,
            OrangeRefGenomeVersion orangeRefGenomeVersion) throws IOException
    {
        ClinicalTranscriptsModel clinicalTranscriptsModel = clinicalTranscriptsTsv != null ?
                ClinicalTranscriptFile.buildFromTsv(orangeRefGenomeVersion, clinicalTranscriptsTsv) : null;
        return new FindingConfig(clinicalTranscriptsModel);
    }

    public FindingConfig(@Nullable final ClinicalTranscriptsModel clinicalTranscriptsModel)
    {
        this.clinicalTranscriptsModel = clinicalTranscriptsModel;
    }

    @Nullable
    public String findCanonicalTranscriptForGene(String gene)
    {
        return clinicalTranscriptsModel != null ? clinicalTranscriptsModel.findCanonicalTranscriptForGene(gene) : null;
    }
}
