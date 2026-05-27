package com.hartwig.hmftools.finding;

import java.io.IOException;
import java.nio.file.Path;

import com.hartwig.hmftools.datamodel.orange.OrangeRefGenomeVersion;
import com.hartwig.hmftools.finding.clinicaltranscript.ClinicalTranscriptFile;
import com.hartwig.hmftools.finding.clinicaltranscript.ClinicalTranscriptsModel;

import org.jetbrains.annotations.Nullable;

record FindingConfig(@Nullable ClinicalTranscriptsModel clinicalTranscriptsModel)
{
    public static FindingConfig createFindingConfig(@Nullable Path clinicalTranscriptsTsv,
            OrangeRefGenomeVersion orangeRefGenomeVersion) throws IOException
    {
        ClinicalTranscriptsModel clinicalTranscriptsModel = clinicalTranscriptsTsv != null ?
                ClinicalTranscriptFile.buildFromTsv(orangeRefGenomeVersion, clinicalTranscriptsTsv) : null;
        return new FindingConfig(clinicalTranscriptsModel);
    }

    @Nullable
    public String findCanonicalTranscriptForGene(String gene)
    {
        return clinicalTranscriptsModel != null ? clinicalTranscriptsModel.findCanonicalTranscriptForGene(gene) : null;
    }
}
