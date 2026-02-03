package com.hartwig.hmftools.common.hla;

import static com.hartwig.hmftools.common.utils.file.FileDelimiters.ITEM_DELIM;

import java.io.File;
import java.io.IOException;
import java.util.List;
import java.util.StringJoiner;

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

    private static final Logger LOGGER = LogManager.getLogger(LilacSummaryData.class);

    public static LilacSummaryData read(final String basePath, final String sampleId) throws IOException
    {
        return load(LilacQcData.generateFilename(basePath, sampleId), LilacAllele.generateFilename(basePath, sampleId));
    }

    public static LilacSummaryData load(final String lilacQcFile, final String lilacResultFile) throws IOException
    {
        LOGGER.info("Loading LILAC data from {}", new File(lilacQcFile).getParent());

        List<LilacQcData> qcData = LilacQcData.read(lilacQcFile);
        StringJoiner qcValues = new StringJoiner(ITEM_DELIM);

        for(LilacQcData geneDataEntry : qcData)
        {
            LOGGER.debug(" Read QC status '{}' for genes '{}' from {}", geneDataEntry.status(), geneDataEntry.genes(), lilacQcFile);
            qcValues.add(geneDataEntry.status());
        }

        List<LilacAllele> alleles = LilacAllele.read(lilacResultFile);
        LOGGER.info(" Read {} LILAC alleles from {}", alleles.size(), lilacResultFile);

        return ImmutableLilacSummaryData.builder()
                .qc(qcValues.toString())
                .alleles(alleles)
                .build();
    }
}
