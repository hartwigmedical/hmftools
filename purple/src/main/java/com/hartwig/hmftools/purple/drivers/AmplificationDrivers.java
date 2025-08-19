package com.hartwig.hmftools.purple.drivers;

import static com.hartwig.hmftools.common.driver.DriverCatalogFactory.createCopyNumberDriver;
import static com.hartwig.hmftools.common.purple.PurpleCommon.DEFAULT_DRIVER_AMPLIFICATION_PLOIDY_RATIO;

import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.stream.Collectors;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.driver.DriverCatalog;
import com.hartwig.hmftools.common.driver.DriverType;
import com.hartwig.hmftools.common.driver.LikelihoodMethod;
import com.hartwig.hmftools.common.driver.panel.DriverGene;
import com.hartwig.hmftools.common.driver.panel.DriverGenePanel;
import com.hartwig.hmftools.common.genome.chromosome.HumanChromosome;
import com.hartwig.hmftools.common.purple.Gender;
import com.hartwig.hmftools.common.purple.PurpleQCStatus;
import com.hartwig.hmftools.common.purple.GeneCopyNumber;
import com.hartwig.hmftools.common.purple.ReportableStatus;

public final class AmplificationDrivers
{
    private static final List<String> TARGET_REGIONS_PARTIAL_AMP_GENES = Lists.newArrayList(
            "BRAF", "EGFR", "CTNNB1", "CBL", "MET", "ALK", "PDGFRA");

    public static List<DriverCatalog> findAmplifications(
            final Set<PurpleQCStatus> qcStatus, final Gender gender, final DriverGenePanel panel,
            final double ploidy, final List<GeneCopyNumber> geneCopyNumbers, boolean isTargetRegions)
    {
        List<DriverCatalog> result = Lists.newArrayList();

        Map<String,DriverGene> amplificationDriverGenes = panel.amplificationTargets().stream()
                .filter(x -> x.reportAmplification())
                .collect(Collectors.toMap(DriverGene::gene, x -> x));

        boolean isHighCopyNoise = qcStatus.contains(PurpleQCStatus.WARN_HIGH_COPY_NUMBER_NOISE);

        for(GeneCopyNumber geneCopyNumber : geneCopyNumbers)
        {
            DriverGene driverGene = amplificationDriverGenes.get(geneCopyNumber.geneName());

            if(isHighCopyNoise && !supportedByOneSV(geneCopyNumber))
                continue;

            double ampCopyNumberThreshold = driverGene != null ?
                    ploidy * driverGene.amplificationRatio() : DEFAULT_DRIVER_AMPLIFICATION_PLOIDY_RATIO;

            double geneCopyNummberThreshold = (gender == Gender.MALE) && HumanChromosome._X.matches(geneCopyNumber.chromosome())
                    ? ampCopyNumberThreshold * 0.5 :  ampCopyNumberThreshold;

            if(geneCopyNumber.minCopyNumber() > geneCopyNummberThreshold)
            {
                geneCopyNumber.setDriverType(DriverType.AMP);
                geneCopyNumber.setReportableStatus(driverGene != null ? ReportableStatus.REPORTED : ReportableStatus.CANDIDATE);
                result.add(createAmpDriver(driverGene, geneCopyNumber));
                continue;
            }

            // partial AMPs are only allowed for WGS
            if(!isTargetRegions || TARGET_REGIONS_PARTIAL_AMP_GENES.contains(geneCopyNumber.geneName()))
            {
                if(geneCopyNumber.maxCopyNumber() > geneCopyNummberThreshold)
                {
                    geneCopyNumber.setDriverType(DriverType.PARTIAL_AMP);
                    geneCopyNumber.setReportableStatus(driverGene != null ? ReportableStatus.REPORTED : ReportableStatus.CANDIDATE);
                    result.add(createPartialAmpDriver(driverGene, geneCopyNumber));
                }
            }
        }

        return result;
    }

    private static boolean supportedByOneSV(final GeneCopyNumber geneCopyNumber)
    {
        return geneCopyNumber.MinRegionStartSupport.isSV() || geneCopyNumber.MinRegionEndSupport.isSV();
    }

    private static DriverCatalog createAmpDriver(final DriverGene driverGene, final GeneCopyNumber geneCopyNumber)
    {
        return createCopyNumberDriver(driverGene.likelihoodType(), DriverType.AMP, LikelihoodMethod.AMP, false, geneCopyNumber);
    }

    private static DriverCatalog createPartialAmpDriver(final DriverGene driverGene, final GeneCopyNumber geneCopyNumber)
    {
        return createCopyNumberDriver(driverGene.likelihoodType(), DriverType.PARTIAL_AMP, LikelihoodMethod.AMP, false, geneCopyNumber);
    }
}
