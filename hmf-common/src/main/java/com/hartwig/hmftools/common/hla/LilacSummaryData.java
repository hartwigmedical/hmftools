package com.hartwig.hmftools.common.hla;

import java.io.File;
import java.io.IOException;
import java.util.List;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.immutables.value.Value;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

@Value.Immutable
@Value.Style(passAnnotations = { NotNull.class, Nullable.class })
public abstract class LilacSummaryData
{
    @NotNull
    public abstract List<LilacQcData> qc();

    @NotNull
    public abstract List<LilacAllele> alleles();

    public int somaticVariantCount()
    {
        return (int)alleles().stream().mapToDouble(x -> x.somaticVariantCount()).sum();
    }

    private static final Logger LOGGER = LogManager.getLogger(LilacSummaryData.class);

    public static LilacSummaryData read(final String basePath, final String sampleId) throws IOException
    {
        List<LilacQcData> qcData = LilacQcData.read(LilacQcData.generateFilename(basePath, sampleId));
        List<LilacAllele> alleles = LilacAllele.read(LilacAllele.generateFilename(basePath, sampleId));

        return ImmutableLilacSummaryData.builder()
                .qc(qcData)
                .alleles(alleles)
                .build();
    }

    public static LilacSummaryData load(final String lilacQcFile, final String lilacResultFile) throws IOException
    {
        LOGGER.info("Loading LILAC data from {}", new File(lilacQcFile).getParent());

        List<LilacQcData> qcData = LilacQcData.read(lilacQcFile);

        for(LilacQcData geneDataEntry : qcData)
            LOGGER.info(" Read QC status '{}' for genes '{}' from {}", geneDataEntry.status(), geneDataEntry.genes(), lilacQcFile);


        List<LilacAllele> alleles = LilacAllele.read(lilacResultFile);
        LOGGER.info(" Read {} LILAC alleles from {}", alleles.size(), lilacResultFile);

        return ImmutableLilacSummaryData.builder()
                .qc(qcData)
                .alleles(alleles)
                .build();
    }
}
