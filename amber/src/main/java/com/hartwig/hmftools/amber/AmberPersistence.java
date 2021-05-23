package com.hartwig.hmftools.amber;

import static com.hartwig.hmftools.amber.AmberConfig.AMB_LOGGER;

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

import org.jetbrains.annotations.NotNull;

class AmberPersistence
{
    private final AmberConfig mConfig;

    public AmberPersistence(final AmberConfig config)
    {
        mConfig = config;
    }

    void persistVersionInfo(@NotNull final VersionInfo versionInfo) throws IOException
    {
        versionInfo.write(mConfig.OutputDir);
    }

    void persistBAF(@NotNull final List<AmberBAF> result) throws IOException, InterruptedException
    {
        final String filename = AmberBAFFile.generateAmberFilenameForWriting(mConfig.OutputDir, mConfig.TumorId);
        AmberBAFFile.write(filename, result);

        AMB_LOGGER.info("Applying pcf segmentation");
        new BAFSegmentation(mConfig.OutputDir).applySegmentation(mConfig.TumorId);
    }

    void persistBafVcf(@NotNull final List<TumorBAF> tumorBAFList, final AmberHetNormalEvidence amberHetNormalEvidence)
    {
        final String outputVcf = mConfig.OutputDir + File.separator + mConfig.TumorId + ".amber.baf.vcf.gz";
        AMB_LOGGER.info("Writing {} BAF records to {}", tumorBAFList.size(), outputVcf);
        new AmberVCF(mConfig).writeBAF(outputVcf, tumorBAFList, amberHetNormalEvidence);
    }

    void persistQC(@NotNull final List<AmberBAF> result, @NotNull final List<TumorContamination> contaminationRecords) throws IOException
    {
        final double contamination = new TumorContaminationModel().contamination(contaminationRecords);
        final AmberQC qcStats = AmberQCFactory.create(contamination, result);
        final String qcFilename = AmberQCFile.generateFilename(mConfig.OutputDir, mConfig.TumorId);
        AmberQCFile.write(qcFilename, qcStats);
    }

    void persistContamination(@NotNull final List<TumorContamination> contaminationList) throws IOException
    {
        Collections.sort(contaminationList);

        final String outputVcf = mConfig.OutputDir + File.separator + mConfig.TumorId + ".amber.contamination.vcf.gz";
        AMB_LOGGER.info("Writing {} contamination records to {}", contaminationList.size(), outputVcf);
        new AmberVCF(mConfig).writeContamination(outputVcf, contaminationList);

        final String filename = TumorContaminationFile.generateContaminationFilename(mConfig.OutputDir, mConfig.TumorId);
        TumorContaminationFile.write(filename, contaminationList);
    }

    void persistSnpCheck(@NotNull final ListMultimap<Chromosome, BaseDepth> baseDepths)
    {
        if(baseDepths.size() > 0)
        {
            final String outputVcf = mConfig.OutputDir + File.separator + mConfig.primaryReference() + ".amber.snp.vcf.gz";
            AMB_LOGGER.info("Writing {} germline snp records to {}", baseDepths.size(), outputVcf);
            new AmberVCF(mConfig).writeSNPCheck(outputVcf, Lists.newArrayList(baseDepths.values()));
        }
    }
}
