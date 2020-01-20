package com.hartwig.hmftools.protect.common;

import java.io.IOException;
import java.util.List;

import com.hartwig.hmftools.common.purple.gene.GeneCopyNumber;
import com.hartwig.hmftools.common.purple.gene.GeneCopyNumberFile;
import com.hartwig.hmftools.common.purple.purity.FittedPurityFile;
import com.hartwig.hmftools.common.purple.purity.PurityContext;
import com.hartwig.hmftools.common.variant.SomaticVariantFactory;
import com.hartwig.hmftools.common.variant.Variant;
import com.hartwig.hmftools.common.variant.structural.annotation.ReportableGeneFusion;
import com.hartwig.hmftools.common.variant.structural.annotation.ReportableGeneFusionFile;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.jetbrains.annotations.NotNull;

public class GenomicData {

    private static final Logger LOGGER = LogManager.getLogger(GenomicData.class);

    private GenomicData() {

    }

    public static double extractPloidy(@NotNull String purplePurityTsv) throws IOException {
        LOGGER.info("Reading purple purity from {}", purplePurityTsv);
        PurityContext purityContext = FittedPurityFile.read(purplePurityTsv);
        double ploidy = purityContext.bestFit().ploidy();
        LOGGER.info(" Sample ploidy: {}", ploidy);
        return ploidy;
    }

    @NotNull
    public static List<? extends Variant> readPassSomaticVariants(@NotNull String sampleId, @NotNull String somaticVariantVcf)
            throws IOException {
        LOGGER.info("Reading somatic variants from {}", somaticVariantVcf);
        List<? extends Variant> passSomaticVariants = SomaticVariantFactory.passOnlyInstance().fromVCFFile(sampleId, somaticVariantVcf);
        LOGGER.info(" Loaded {} PASS somatic variants", passSomaticVariants.size());
        return passSomaticVariants;
    }

    @NotNull
    public static List<GeneCopyNumber> readGeneCopyNumbers(@NotNull String purpleGeneCnvTsv) throws IOException {
        LOGGER.info("Reading gene copy numbers from {}", purpleGeneCnvTsv);
        List<GeneCopyNumber> geneCopyNumbers = GeneCopyNumberFile.read(purpleGeneCnvTsv);
        LOGGER.info(" Loaded {} gene copy numbers", geneCopyNumbers.size());
        return geneCopyNumbers;
    }

    @NotNull
    public static List<ReportableGeneFusion> readGeneFusions(@NotNull String linxFusionTsv) throws IOException {
        LOGGER.info("Reading gene fusions from {}", linxFusionTsv);
        List<ReportableGeneFusion> fusions = ReportableGeneFusionFile.read(linxFusionTsv);
        LOGGER.info(" Loaded {} fusions", fusions.size());
        return fusions;
    }
}
