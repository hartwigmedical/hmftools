package com.hartwig.hmftools.purple.drivers;

import static com.hartwig.hmftools.common.driver.DriverCatalogFactory.createCopyNumberDriver;
import static com.hartwig.hmftools.common.driver.DriverType.DEL;
import static com.hartwig.hmftools.common.driver.DriverType.HET_DEL;
import static com.hartwig.hmftools.common.driver.DriverType.UNKNOWN;

import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.stream.Collectors;

import com.google.common.collect.Lists;
import com.google.common.collect.Sets;
import com.hartwig.hmftools.common.driver.DriverCatalog;
import com.hartwig.hmftools.common.driver.DriverType;
import com.hartwig.hmftools.common.driver.LikelihoodMethod;
import com.hartwig.hmftools.common.driver.panel.DriverGene;
import com.hartwig.hmftools.common.driver.panel.DriverGenePanel;
import com.hartwig.hmftools.common.genome.chromosome.HumanChromosome;
import com.hartwig.hmftools.common.purple.GeneCopyNumber;
import com.hartwig.hmftools.common.purple.PurpleQCStatus;
import com.hartwig.hmftools.common.purple.ReportableStatus;
import com.hartwig.hmftools.common.purple.SegmentSupport;
import com.hartwig.hmftools.common.utils.Doubles;

public class DeletionDrivers
{
    public static final double MAX_COPY_NUMBER_DEL = 0.5;

    private static final int SHORT_DEL_LENGTH = 10_000_000;
    private static final double SINGLE_DEPTH_GC_MIN = 0.35;
    private static final double SINGLE_DEPTH_GC_MAX = 0.6;

    private static final Set<SegmentSupport> MERE = Sets.newHashSet(SegmentSupport.CENTROMERE, SegmentSupport.TELOMERE);

    public static int deletedGenes(final List<GeneCopyNumber> geneCopyNumbers)
    {
        return (int) geneCopyNumbers.stream()
                .filter(x -> !HumanChromosome.fromString(x.chromosome()).equals(HumanChromosome._Y)
                        && Doubles.lessThan(x.minCopyNumber(), MAX_COPY_NUMBER_DEL))
                .count();
    }

    public static List<DriverCatalog> findDeletions(
            final Set<PurpleQCStatus> qcStatus, final double ploidy, final DriverGenePanel panel,
            final List<GeneCopyNumber> geneCopyNumbers, boolean isTargetRegions)
    {
        List<DriverCatalog> drivers = Lists.newArrayList();

        boolean checkQcStatus = qcStatus.contains(PurpleQCStatus.WARN_DELETED_GENES)
                || qcStatus.contains(PurpleQCStatus.WARN_HIGH_COPY_NUMBER_NOISE);

        Map<String,DriverGene> deletionGenes = panel.deletionTargets().stream()
                .filter(x -> x.reportDeletion() || x.reportHetDeletion())
                .collect(Collectors.toMap(DriverGene::gene, x -> x));

        for(GeneCopyNumber geneCopyNumber : geneCopyNumbers)
        {
            DriverGene driverGene = deletionGenes.get(geneCopyNumber.geneName());

            if(driverGene == null)
                continue;

            DriverType driverType = DriverType.UNKNOWN;

            if(driverGene.reportDeletion() && geneCopyNumber.minCopyNumber() < MAX_COPY_NUMBER_DEL)
            {
                driverType = DEL;
            }
            else if(driverGene.reportHetDeletion() && ploidy > 0 && HumanChromosome.fromString(geneCopyNumber.Chromosome).isAutosome())
            {
                double adjustedMinCopyNumber = geneCopyNumber.minCopyNumber() / ploidy;

                if(adjustedMinCopyNumber < driverGene.hetDeletionThreshold())
                    driverType = HET_DEL;
            }

            if(driverType == UNKNOWN)
                continue;

            boolean hasFullSvSupport = supportedByTwoSVs(geneCopyNumber);
            boolean isShort = geneCopyNumber.minRegionBases() < SHORT_DEL_LENGTH;

            if(isTargetRegions)
            {
                if(!hasFullSvSupport)
                {
                    if(checkQcStatus)
                        continue;

                    if(geneCopyNumber.DepthWindowCount == 1)
                    {
                        if(geneCopyNumber.GcContent < SINGLE_DEPTH_GC_MIN || geneCopyNumber.GcContent > SINGLE_DEPTH_GC_MAX)
                            continue;
                    }
                }
            }
            else
            {
                if(checkQcStatus)
                {
                    if(!hasFullSvSupport && !(isShort && supportedByOneSVAndMere(geneCopyNumber)))
                        continue;
                }
            }

            geneCopyNumber.setDriverType(driverType);
            geneCopyNumber.setReportableStatus(ReportableStatus.REPORTED); // candidates not yet supported
            drivers.add(createDelDriver(driverGene, driverType, geneCopyNumber));
        }

        return drivers;
    }

    private static boolean supportedByTwoSVs(final GeneCopyNumber geneCopyNumber)
    {
        return geneCopyNumber.MinRegionStartSupport.isSV() && geneCopyNumber.MinRegionEndSupport.isSV();
    }

    private static boolean supportedByOneSVAndMere(final GeneCopyNumber geneCopyNumber)
    {
        if(MERE.contains(geneCopyNumber.MinRegionStartSupport) && geneCopyNumber.MinRegionEndSupport.isSV())
        {
            return true;
        }

        if(MERE.contains(geneCopyNumber.MinRegionEndSupport) && geneCopyNumber.MinRegionStartSupport.isSV())
        {
            return true;
        }

        return false;
    }

    private static DriverCatalog createDelDriver(
            final DriverGene driverGene, final DriverType driverType, final GeneCopyNumber geneCopyNumber)
    {
        boolean biallelic = driverType == DEL;
        double likelihood = driverType == DEL ? 1 : 0;
        return createCopyNumberDriver(driverGene.likelihoodType(), driverType, LikelihoodMethod.DEL, biallelic, likelihood, geneCopyNumber);
    }
}
