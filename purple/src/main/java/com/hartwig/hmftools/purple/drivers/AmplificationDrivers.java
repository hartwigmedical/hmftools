package com.hartwig.hmftools.purple.drivers;

import static com.hartwig.hmftools.common.driver.DriverType.UNKNOWN;
import static com.hartwig.hmftools.purple.drivers.DeletionDrivers.createCopyNumberDriver;

import java.util.List;
import java.util.Set;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.driver.DriverCatalog;
import com.hartwig.hmftools.common.driver.DriverCategory;
import com.hartwig.hmftools.common.driver.DriverType;
import com.hartwig.hmftools.common.driver.LikelihoodMethod;
import com.hartwig.hmftools.common.driver.panel.DriverGene;
import com.hartwig.hmftools.purple.DriverGeneResource;
import com.hartwig.hmftools.common.genome.chromosome.HumanChromosome;
import com.hartwig.hmftools.common.purple.Gender;
import com.hartwig.hmftools.common.purple.PurpleQCStatus;
import com.hartwig.hmftools.common.purple.GeneCopyNumber;
import com.hartwig.hmftools.common.purple.ReportedStatus;

public final class AmplificationDrivers
{
    private static final List<String> TARGET_REGIONS_PARTIAL_AMP_GENES = Lists.newArrayList(
            "BRAF", "EGFR", "CTNNB1", "CBL", "MET", "ALK", "PDGFRA");

    public static List<DriverCatalog> findAmplifications(
            final Set<PurpleQCStatus> qcStatus, final Gender gender, final DriverGeneResource panel,
            final double ploidy, final List<GeneCopyNumber> geneCopyNumbers, boolean isTargetRegions)
    {
        List<DriverCatalog> result = Lists.newArrayList();

        boolean isHighCopyNoise = qcStatus.contains(PurpleQCStatus.WARN_HIGH_COPY_NUMBER_NOISE);

        for(GeneCopyNumber geneCopyNumber : geneCopyNumbers)
        {
            DriverGene driverGene = panel.DriverGeneMap.get(geneCopyNumber.geneName());

            if(driverGene == null)
                continue;

            if(isHighCopyNoise && !supportedByOneSV(geneCopyNumber))
                continue;

            double ampCopyNumberThreshold = ploidy * driverGene.amplificationRatio();

            double geneCopyNummberThreshold = (gender == Gender.MALE) && HumanChromosome._X.matches(geneCopyNumber.chromosome())
                    ? ampCopyNumberThreshold * 0.5 :  ampCopyNumberThreshold;

            DriverType driverType = UNKNOWN;

            if(geneCopyNumber.minCopyNumber() > geneCopyNummberThreshold)
            {
                driverType = DriverType.AMP;
            }
            else if(!isTargetRegions || TARGET_REGIONS_PARTIAL_AMP_GENES.contains(geneCopyNumber.geneName()))
            {
                if(geneCopyNumber.maxCopyNumber() > geneCopyNummberThreshold)
                {
                    driverType = DriverType.PARTIAL_AMP;
                }
            }

            ReportedStatus reportedStatus = driverType != UNKNOWN ? ReportedStatus.REPORTED : ReportedStatus.NONE;

            geneCopyNumber.setDriverType(driverType);
            geneCopyNumber.setReportedStatus(reportedStatus);

            DriverCatalog driverCatalog = createCopyNumberDriver(
                    DriverCategory.ONCO, driverType, LikelihoodMethod.AMP, false, 1, geneCopyNumber, reportedStatus);
            result.add(driverCatalog);
        }

        return result;
    }

    private static boolean supportedByOneSV(final GeneCopyNumber geneCopyNumber)
    {
        return geneCopyNumber.MinRegionStartSupport.isSV() || geneCopyNumber.MinRegionEndSupport.isSV();
    }
}
