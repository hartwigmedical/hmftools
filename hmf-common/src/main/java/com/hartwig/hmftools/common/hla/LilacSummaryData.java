package com.hartwig.hmftools.common.hla;

import java.io.File;
import java.io.IOException;
import java.util.List;
import java.util.Map;

import com.google.common.collect.Maps;

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
    public abstract Map<String, String> qc__();

    @NotNull
    public abstract List<LilacAllele> alleles__();

    // TODO: by genes?
//    public int somaticVariantCount_()
//    {
//        return (int)alleles_().stream().mapToDouble(x -> x.somaticVariantCount()).sum();
//    }

    private static final Logger LOGGER = LogManager.getLogger(LilacSummaryData.class);

    public static LilacSummaryData read_(final String basePath, final String sampleId) throws IOException
    {
        List<LilacQcData> qcData_ = LilacQcData.read_(LilacQcData.generateFilename(basePath, sampleId));
        List<LilacAllele> alleles = LilacAllele.read(LilacAllele.generateFilename(basePath, sampleId));

        // TODO: common code?
        Map<String, String> qcStatusByGenes = Maps.newHashMap();
        for(LilacQcData qcDataEntry : qcData_)
            qcStatusByGenes.put(qcDataEntry.genes(), qcDataEntry.status());

        return ImmutableLilacSummaryData.builder()
                .qc__(qcStatusByGenes)
                .alleles__(alleles)
                .build();
    }

    public static LilacSummaryData load_(final String lilacQcFile, final String lilacResultFile) throws IOException
    {
        LOGGER.info("Loading LILAC data from {}", new File(lilacQcFile).getParent());

        List<LilacQcData> qcData_ = LilacQcData.read_(lilacQcFile);

        // TODO: common code?
        Map<String, String> qcStatusByGenes = Maps.newHashMap();
        for(LilacQcData qcDataEntry : qcData_)
            qcStatusByGenes.put(qcDataEntry.genes(), qcDataEntry.status());

        LOGGER.info(" Read QC status by genes '{}' from {}", qcStatusByGenes, lilacQcFile);

        List<LilacAllele> alleles = LilacAllele.read(lilacResultFile);

        LOGGER.info(" Read {} LILAC alleles from {}", alleles.size(), lilacResultFile);

        return ImmutableLilacSummaryData.builder()
                .qc__(qcStatusByGenes)
                .alleles__(alleles)
                .build();
    }
}
