package com.hartwig.hmftools.amber;

import static com.hartwig.hmftools.amber.AmberConfig.AMB_LOGGER;

import java.io.IOException;
import java.util.Collections;
import java.util.List;

import com.google.common.collect.ListMultimap;
import com.hartwig.hmftools.common.amber.AmberBAF;
import com.hartwig.hmftools.common.amber.AmberBAFFile;
import com.hartwig.hmftools.common.amber.AmberQC;
import com.hartwig.hmftools.common.amber.AmberQCFile;
import com.hartwig.hmftools.common.amber.ImmutableAmberQC;
import com.hartwig.hmftools.common.genome.chromosome.Chromosome;
import com.hartwig.hmftools.common.utils.version.VersionInfo;

import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

public class ResultsWriter
{
    private final AmberConfig mConfig;

    public ResultsWriter(final AmberConfig config)
    {
        mConfig = config;
    }

    void persistVersionInfo(final VersionInfo versionInfo) throws IOException
    {
        versionInfo.write(mConfig.OutputDir);
    }

    void persistBAF(final List<AmberBAF> result) throws IOException, InterruptedException
    {
        final String filename = AmberBAFFile.generateAmberFilenameForWriting(mConfig.OutputDir, mConfig.getSampleId());
        AmberBAFFile.write(filename, result);

        if(mConfig.TumorId != null && !mConfig.SkipBafSegmentation)
        {
            AMB_LOGGER.info("applying pcf segmentation");
            new BAFSegmentation(mConfig.OutputDir).applySegmentation(mConfig.TumorId, filename);
        }
    }

    void persistQC(final List<TumorContamination> contaminationRecords,
            double consanguinityProportion, @Nullable Chromosome uniparentalDisomy) throws IOException
    {
        final double contamination = new TumorContaminationModel().contamination(contaminationRecords);

        AmberQC qcStats = ImmutableAmberQC.builder()
                .contamination(contamination)
                .consanguinityProportion(consanguinityProportion)
                .uniparentalDisomy(uniparentalDisomy != null ? uniparentalDisomy.toString() : null).build();

        final String qcFilename = AmberQCFile.generateFilename(mConfig.OutputDir, mConfig.getSampleId());
        AmberQCFile.write(qcFilename, qcStats);
    }

    void persistContamination(final List<TumorContamination> contaminationList) throws IOException
    {
        Collections.sort(contaminationList);

        final String outputVcf = mConfig.OutputDir + mConfig.TumorId + ".amber.contamination.vcf.gz";
        AMB_LOGGER.info("writing {} contamination records to {}", contaminationList.size(), outputVcf);
        new VCFWriter(mConfig).writeContamination(outputVcf, contaminationList);

        final String filename = TumorContaminationFile.generateContaminationFilename(mConfig.OutputDir, mConfig.TumorId);
        TumorContaminationFile.write(filename, contaminationList);
    }

    void persistSnpCheck(final ListMultimap<Chromosome, PositionEvidence> baseDepths)
    {
        if(baseDepths.size() > 0)
        {
            final String outputVcf = mConfig.OutputDir + mConfig.primaryReference() + ".amber.snp.vcf.gz";
            AMB_LOGGER.info("writing {} germline snp records to {}", baseDepths.size(), outputVcf);
            VCFWriter.writeBaseDepths(outputVcf, baseDepths.values(), mConfig.primaryReference());
        }
    }

    void persistHomozygousRegions(final List<RegionOfHomozygosity> regionOfHomozygosities) throws IOException
    {
        final String filename = RegionOfHomozygosityFile.generateFilename(mConfig.OutputDir, mConfig.primaryReference());
        RegionOfHomozygosityFile.write(filename, regionOfHomozygosities);
    }
}