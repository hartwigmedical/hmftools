package com.hartwig.hmftools.amber;

import java.io.File;
import java.io.IOException;
import java.util.Collections;
import java.util.List;

import com.hartwig.hmftools.common.amber.AmberBAF;
import com.hartwig.hmftools.common.amber.AmberBAFFile;
import com.hartwig.hmftools.common.amber.AmberVCF;
import com.hartwig.hmftools.common.amber.TumorBAF;
import com.hartwig.hmftools.common.amber.TumorContamination;
import com.hartwig.hmftools.common.amber.TumorContaminationFile;
import com.hartwig.hmftools.common.amber.TumorContaminationModel;
import com.hartwig.hmftools.common.amber.qc.AmberQC;
import com.hartwig.hmftools.common.amber.qc.AmberQCFactory;
import com.hartwig.hmftools.common.amber.qc.AmberQCFile;
import com.hartwig.hmftools.common.version.VersionInfo;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.jetbrains.annotations.NotNull;

class AmberPersistence {

    private static final Logger LOGGER = LogManager.getLogger(AmberPersistence.class);

    private final AmberConfig config;

    AmberPersistence(final AmberConfig config) {
        this.config = config;
    }

    void persistVersionInfo() throws IOException {

        final VersionInfo versionInfo = new VersionInfo("amber.version");
        versionInfo.write(config.outputDirectory());
    }

    void persistAmberBAF(@NotNull final List<AmberBAF> result) throws IOException, InterruptedException {
        LOGGER.info("Writing {} BAF records to {}", result.size(), config.outputDirectory());
        final String filename = AmberBAFFile.generateAmberFilename(config.outputDirectory(), config.tumor());
        AmberBAFFile.write(filename, result);

        LOGGER.info("Applying pcf segmentation");
        new BAFSegmentation(config.outputDirectory()).applySegmentation(config.tumor());
    }

    void persistTumorBAF(@NotNull final List<TumorBAF> tumorBAFList) {
        final String outputVcf = config.outputDirectory() + File.separator + config.tumor() + ".amber.vcf.gz";
        new AmberVCF(config.normal(), config.tumor()).write(outputVcf, tumorBAFList);
    }

    void persisQC(@NotNull final List<AmberBAF> result, @NotNull final List<TumorContamination> contaminationRecords) throws IOException {
        final double contamination = new TumorContaminationModel().contamination(contaminationRecords);
        final AmberQC qcStats = AmberQCFactory.create(contamination, result);
        final String qcFilename = AmberQCFile.generateFilename(config.outputDirectory(), config.tumor());
        AmberQCFile.write(qcFilename, qcStats);
    }

    void persistContamination(@NotNull final List<TumorContamination> contaminationList) throws IOException {
        Collections.sort(contaminationList);

        LOGGER.info("Writing {} contamination records to {}", contaminationList.size(), config.outputDirectory());
        final String outputVcf = config.outputDirectory() + File.separator + config.tumor() + ".amber.contamination.vcf.gz";
        new AmberVCF(config.normal(), config.tumor()).writeContamination(outputVcf, contaminationList);

        final String filename = TumorContaminationFile.generateContaminationFilename(config.outputDirectory(), config.tumor());
        TumorContaminationFile.write(filename, contaminationList);
    }

}
