package com.hartwig.hmftools.common.hla;

import static com.hartwig.hmftools.common.utils.file.FileDelimiters.ITEM_DELIM;

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

    public static LilacSummaryData load(final String qcFile, final String resultsFile) throws IOException
    {
        List<LilacQcData> qcData = LilacQcData.read(qcFile);
        StringJoiner qcValues = new StringJoiner(ITEM_DELIM);

        for(LilacQcData geneDataEntry : qcData)
        {
            LOGGER.debug(" Read QC status '{}' for genes '{}' from {}", geneDataEntry.status(), geneDataEntry.genes(), qcFile);
            qcValues.add(geneDataEntry.status());
        }

        List<LilacAllele> alleles = LilacAllele.read(resultsFile);
        LOGGER.info(" Read {} Lilac alleles from {}", alleles.size(), resultsFile);

        return ImmutableLilacSummaryData.builder()
                .qc(qcValues.toString())
                .alleles(alleles)
                .build();
    }
}
