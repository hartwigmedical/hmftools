package com.hartwig.hmftools.amber;

import java.io.File;
import java.io.IOException;
import java.util.Collections;
import java.util.List;

import com.google.common.collect.ListMultimap;
import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.amber.AmberBAF;
import com.hartwig.hmftools.common.amber.AmberBAFFile;
import com.hartwig.hmftools.common.amber.BaseDepth;
import com.hartwig.hmftools.common.amber.TumorBAF;
import com.hartwig.hmftools.common.amber.TumorContamination;
import com.hartwig.hmftools.common.amber.TumorContaminationFile;
import com.hartwig.hmftools.common.amber.TumorContaminationModel;
import com.hartwig.hmftools.common.amber.qc.AmberQC;
import com.hartwig.hmftools.common.amber.qc.AmberQCFactory;
import com.hartwig.hmftools.common.amber.qc.AmberQCFile;
import com.hartwig.hmftools.common.genome.chromosome.Chromosome;
import com.hartwig.hmftools.common.utils.version.VersionInfo;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.jetbrains.annotations.NotNull;

class AmberPersistence {

    private static final Logger LOGGER = LogManager.getLogger(AmberPersistence.class);

    private final AmberConfig config;

    AmberPersistence(final AmberConfig config) {
        this.config = config;
    }

    void persistVersionInfo(@NotNull final VersionInfo versionInfo) throws IOException {
        versionInfo.write(config.outputDirectory());
    }

    void persistBAF(@NotNull final List<AmberBAF> result) throws IOException, InterruptedException {
        final String filename = AmberBAFFile.generateAmberFilenameForWriting(config.outputDirectory(), config.tumor());
        AmberBAFFile.write(filename, result);

        LOGGER.info("Applying pcf segmentation");
        new BAFSegmentation(config.outputDirectory()).applySegmentation(config.tumor());
    }


    void persistBafVcf(@NotNull final List<TumorBAF> tumorBAFList, final AmberHetNormalEvidence amberHetNormalEvidence) {
        final String outputVcf = config.outputDirectory() + File.separator + config.tumor() + ".amber.baf.vcf.gz";
        LOGGER.info("Writing {} BAF records to {}", tumorBAFList.size(), outputVcf);
        new AmberVCF(config).writeBAF(outputVcf, tumorBAFList, amberHetNormalEvidence);
    }

    void persistQC(@NotNull final List<AmberBAF> result, @NotNull final List<TumorContamination> contaminationRecords) throws IOException {
        final double contamination = new TumorContaminationModel().contamination(contaminationRecords);
        final AmberQC qcStats = AmberQCFactory.create(contamination, result);
        final String qcFilename = AmberQCFile.generateFilename(config.outputDirectory(), config.tumor());
        AmberQCFile.write(qcFilename, qcStats);
    }

    void persistContamination(@NotNull final List<TumorContamination> contaminationList) throws IOException {
        Collections.sort(contaminationList);

        final String outputVcf = config.outputDirectory() + File.separator + config.tumor() + ".amber.contamination.vcf.gz";
        LOGGER.info("Writing {} contamination records to {}", contaminationList.size(), outputVcf);
        new AmberVCF(config).writeContamination(outputVcf, contaminationList);

        final String filename = TumorContaminationFile.generateContaminationFilename(config.outputDirectory(), config.tumor());
        TumorContaminationFile.write(filename, contaminationList);
    }

    void persistSnpCheck(@NotNull final ListMultimap<Chromosome, BaseDepth> baseDepths) {
        if (baseDepths.size() > 0) {
            final String outputVcf = config.outputDirectory() + File.separator + config.reference().get(0) + ".amber.snp.vcf.gz";
            LOGGER.info("Writing {} germline snp records to {}", baseDepths.size(), outputVcf);
            new AmberVCF(config).writeSNPCheck(outputVcf, Lists.newArrayList(baseDepths.values()));
        }
    }
}
