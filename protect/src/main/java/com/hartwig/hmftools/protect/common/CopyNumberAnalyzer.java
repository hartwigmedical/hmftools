package com.hartwig.hmftools.protect.common;

import java.io.IOException;
import java.text.DecimalFormat;
import java.util.List;

import com.hartwig.hmftools.common.purple.CheckPurpleQuality;
import com.hartwig.hmftools.common.purple.gene.GeneCopyNumber;
import com.hartwig.hmftools.common.purple.gene.GeneCopyNumberFile;
import com.hartwig.hmftools.common.purple.purity.FittedPurity;
import com.hartwig.hmftools.common.purple.purity.FittedPurityFile;
import com.hartwig.hmftools.common.purple.purity.PurityContext;
import com.hartwig.hmftools.common.purple.qc.PurpleQC;
import com.hartwig.hmftools.common.purple.qc.PurpleQCFile;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.jetbrains.annotations.NotNull;

public class CopyNumberAnalyzer {

    private CopyNumberAnalyzer(){

    }
    private static final Logger LOGGER = LogManager.getLogger(CopyNumberAnalyzer.class);


    @NotNull
    public static CopyNumberAnalysis analyzeCopyNumbers(@NotNull String purplePurityTsv, @NotNull String purpleQCFile,
            @NotNull String purpleGeneCnvTsv) throws IOException {
        PurityContext purityContext = FittedPurityFile.read(purplePurityTsv);
        PurpleQC purpleQC = PurpleQCFile.read(purpleQCFile);

        LOGGER.info("Loaded purple sample data from {}", purplePurityTsv);
        LOGGER.info(" Purple purity: {}", new DecimalFormat("#'%'").format(purityContext.bestFit().purity() * 100));
        LOGGER.info(" Purple average tumor ploidy: {}", purityContext.bestFit().ploidy());
        LOGGER.info(" Purple status: {}", purityContext.status());
        LOGGER.info(" WGD happened: {}", purityContext.wholeGenomeDuplication() ? "yes" : "no");
        LOGGER.info("Loaded purple QC data from {}", purpleQCFile);
        LOGGER.info(" Purple QC status: {}", purpleQC.status());

        List<GeneCopyNumber> exomeGeneCopyNumbers = GeneCopyNumberFile.read(purpleGeneCnvTsv);
        LOGGER.info("Loaded {} gene copy numbers from {}", exomeGeneCopyNumbers.size(), purpleGeneCnvTsv);

        FittedPurity bestFit = purityContext.bestFit();
        List<ReportableGainLoss> reportableGainsAndLosses =
                ExtractReportableGainsAndLosses.toReportableGainsAndLosses(exomeGeneCopyNumbers, bestFit.ploidy());

        return ImmutableCopyNumberAnalysis.builder()
                .purity(bestFit.purity())
                .hasReliablePurity(CheckPurpleQuality.checkHasReliablePurity(purityContext))
                .hasReliableQuality(CheckPurpleQuality.checkHasReliableQuality(purpleQC))
                .ploidy(bestFit.ploidy())
                .exomeGeneCopyNumbers(exomeGeneCopyNumbers)
                .reportableGainsAndLosses(reportableGainsAndLosses)
                .microsatelliteIndelsPerMb(purityContext.microsatelliteIndelsPerMb())
                .microsatelliteStatus(purityContext.microsatelliteStatus())
                .tumorMutationalLoad(purityContext.tumorMutationalLoad())
                .tumorMutationalLoadStatus(purityContext.tumorMutationalLoadStatus())
                .tumorMutationalBurdenPerMb(purityContext.tumorMutationalBurdenPerMb())
                .build();
    }
}
