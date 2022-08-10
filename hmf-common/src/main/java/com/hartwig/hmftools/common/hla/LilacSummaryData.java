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
    public abstract String qc();

    @NotNull
    public abstract List<LilacAllele> alleles();

    public int somaticVariantCount()
    {
        return (int)alleles().stream().mapToDouble(x -> x.somaticVariantCount()).sum();
    }

    private static final Logger LOGGER = LogManager.getLogger(LilacSummaryData.class);

    public static LilacSummaryData read(final String basePath, final String sampleId) throws IOException
    {
        LilacQcData qcData = LilacQcData.read(LilacQcData.generateFilename(basePath, sampleId));
        List<LilacAllele> alleles = LilacAllele.read(LilacAllele.generateFilename(basePath, sampleId));

        return ImmutableLilacSummaryData.builder()
                .qc(qcData.status())
                .alleles(alleles)
                .build();
    }

    public static LilacSummaryData load(final String lilacQcCsv, final String lilacResultCsv) throws IOException
    {
        LOGGER.info("Loading LILAC data from {}", new File(lilacQcCsv).getParent());

        LilacQcData qcData = LilacQcData.read(lilacQcCsv);

        LOGGER.info("Read QC status '{}' from {}", qcData.status(), lilacQcCsv);

        List<LilacAllele> alleles = LilacAllele.read(lilacResultCsv);
        LOGGER.info("Read {} LILAC alleles from {}", alleles.size(), lilacResultCsv);

        return ImmutableLilacSummaryData.builder()
                .qc(qcData.status())
                .alleles(alleles)
                .build();
    }
}
