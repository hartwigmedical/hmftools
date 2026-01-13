package com.hartwig.hmftools.purple.drivers;

import static com.hartwig.hmftools.common.driver.DriverType.UNKNOWN;
import static com.hartwig.hmftools.common.purple.PurpleCommon.DRIVER_AMPLIFICATION_CANDIDATE_PLOIDY_RATIO;
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
                    ? ampCopyNumberThreshold * 0.5 : ampCopyNumberThreshold;

            boolean checkPartials = !isTargetRegions || TARGET_REGIONS_PARTIAL_AMP_GENES.contains(geneCopyNumber.geneName());

            boolean isCandidate = false;
            DriverType driverType = calcDriverType(geneCopyNumber, geneCopyNummberThreshold, checkPartials);

            if(driverType == UNKNOWN)
            {
                double candidateAmpThreshold = geneCopyNummberThreshold * DRIVER_AMPLIFICATION_CANDIDATE_PLOIDY_RATIO;
                driverType = calcDriverType(geneCopyNumber, candidateAmpThreshold, checkPartials);
                isCandidate = driverType != UNKNOWN;
            }

            if(driverType == UNKNOWN)
                continue;

            ReportedStatus reportedStatus = driverGene.reportAmplification() ? ReportedStatus.REPORTED : ReportedStatus.NOT_REPORTED;

            geneCopyNumber.setDriverType(driverType);
            geneCopyNumber.setReportedStatus(reportedStatus);
            double likelihood = isCandidate ? 0 : 1;

            DriverCatalog driverCatalog = createCopyNumberDriver(
                    DriverCategory.ONCO, driverType, LikelihoodMethod.AMP, false, likelihood, geneCopyNumber, reportedStatus);
            result.add(driverCatalog);
        }

        return result;
    }

    private static DriverType calcDriverType(final GeneCopyNumber geneCopyNumber, final double geneCopyNumberThreshold, boolean checkPartials)
    {
        if(geneCopyNumber.minCopyNumber() > geneCopyNumberThreshold)
        {
            return DriverType.AMP;
        }
        else if(checkPartials)
        {
            if(geneCopyNumber.maxCopyNumber() > geneCopyNumberThreshold)
            {
                return DriverType.PARTIAL_AMP;
            }
        }

        return UNKNOWN;
    }

    private static boolean supportedByOneSV(final GeneCopyNumber geneCopyNumber)
    {
        return geneCopyNumber.MinRegionStartSupport.isSV() || geneCopyNumber.MinRegionEndSupport.isSV();
    }
}
