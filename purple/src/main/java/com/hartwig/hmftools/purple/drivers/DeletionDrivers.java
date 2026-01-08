package com.hartwig.hmftools.purple.drivers;

import static com.hartwig.hmftools.common.driver.DriverType.DEL;
import static com.hartwig.hmftools.common.driver.DriverType.HET_DEL;
import static com.hartwig.hmftools.common.driver.DriverType.LOH;
import static com.hartwig.hmftools.common.driver.DriverType.UNKNOWN;
import static com.hartwig.hmftools.common.purple.PurpleCommon.LOH_MINOR_ALLEL_CN;

import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.stream.Collectors;

import com.google.common.collect.Lists;
import com.google.common.collect.Sets;
import com.hartwig.hmftools.common.driver.DriverCatalog;
import com.hartwig.hmftools.common.driver.DriverCategory;
import com.hartwig.hmftools.common.driver.DriverType;
import com.hartwig.hmftools.common.driver.ImmutableDriverCatalog;
import com.hartwig.hmftools.common.driver.LikelihoodMethod;
import com.hartwig.hmftools.common.driver.panel.DriverGene;
import com.hartwig.hmftools.purple.DriverGeneResource;
import com.hartwig.hmftools.common.genome.chromosome.HumanChromosome;
import com.hartwig.hmftools.common.purple.GeneCopyNumber;
import com.hartwig.hmftools.common.purple.PurpleQCStatus;
import com.hartwig.hmftools.common.purple.ReportedStatus;
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
            final Set<PurpleQCStatus> qcStatus, final double ploidy, final DriverGeneResource panel,
            final List<GeneCopyNumber> geneCopyNumbers, boolean isTargetRegions)
    {
        List<DriverCatalog> drivers = Lists.newArrayList();

        boolean checkQcStatus = qcStatus.contains(PurpleQCStatus.WARN_DELETED_GENES)
                || qcStatus.contains(PurpleQCStatus.WARN_HIGH_COPY_NUMBER_NOISE);

        for(GeneCopyNumber geneCopyNumber : geneCopyNumbers)
        {
            DriverGene driverGene = panel.DriverGeneMap.get(geneCopyNumber.geneName());

            if(driverGene == null)
                continue;

            boolean hasFullSvSupport = supportedByTwoSVs(geneCopyNumber);
            boolean isShort = geneCopyNumber.minRegionBases() < SHORT_DEL_LENGTH;

            if(isTargetRegions)
            {
                if(!hasFullSvSupport)
                {
                    if(checkQcStatus)
                        continue;

                    if(geneCopyNumber.DepthWindowCount == 1
                    && (geneCopyNumber.GcContent < SINGLE_DEPTH_GC_MIN || geneCopyNumber.GcContent > SINGLE_DEPTH_GC_MAX))
                        continue;
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

            DriverType driverType = DriverType.UNKNOWN;
            ReportedStatus reportedStatus = ReportedStatus.NONE;

            if(geneCopyNumber.minCopyNumber() < MAX_COPY_NUMBER_DEL)
            {
                driverType = DEL;
                reportedStatus = driverGene.reportDeletion() ? ReportedStatus.REPORTED : ReportedStatus.NOT_REPORTED;
            }
            else if(ploidy > 0 && HumanChromosome.fromString(geneCopyNumber.Chromosome).isAutosome())
            {
                double adjustedMinCopyNumber = geneCopyNumber.minCopyNumber() / ploidy;

                if(adjustedMinCopyNumber < driverGene.hetDeletionThreshold())
                {
                    driverType = HET_DEL;
                    reportedStatus = driverGene.reportHetDeletion() ? ReportedStatus.REPORTED : ReportedStatus.NOT_REPORTED;
                }
            }

            if(driverType == UNKNOWN && geneCopyNumber.MinMinorAlleleCopyNumber < LOH_MINOR_ALLEL_CN)
            {
                driverType = LOH;
                reportedStatus = driverGene.reportLoh() ? ReportedStatus.REPORTED : ReportedStatus.NOT_REPORTED;
            }

            if(driverType == UNKNOWN)
                continue;

            geneCopyNumber.setDriverType(driverType);
            geneCopyNumber.setReportedStatus(reportedStatus);

            boolean biallelic = driverType == DEL;
            double likelihood = driverType == DEL ? 1 : 0;

            DriverCatalog driverCatalog = createCopyNumberDriver(
                    driverGene.likelihoodType(), driverType, LikelihoodMethod.DEL, biallelic, likelihood, geneCopyNumber, reportedStatus);

            drivers.add(driverCatalog);
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

    protected static DriverCatalog createCopyNumberDriver(
            DriverCategory category, DriverType driver, final LikelihoodMethod likelihoodMethod, final boolean biallelic,
            final double likelihood, final GeneCopyNumber geneCopyNumber, final ReportedStatus reportedStatus)
    {
        return ImmutableDriverCatalog.builder()
                .chromosome(geneCopyNumber.chromosome())
                .chromosomeBand(geneCopyNumber.ChromosomeBand)
                .gene(geneCopyNumber.geneName())
                .transcript(geneCopyNumber.TransName)
                .isCanonical(geneCopyNumber.IsCanonical)
                .missense(0)
                .nonsense(0)
                .inframe(0)
                .frameshift(0)
                .splice(0)
                .driverLikelihood(likelihood)
                .driver(driver)
                .likelihoodMethod(likelihoodMethod)
                .category(category)
                .biallelic(biallelic)
                .minCopyNumber(geneCopyNumber.minCopyNumber())
                .maxCopyNumber(geneCopyNumber.maxCopyNumber())
                .reportedStatus(reportedStatus)
                .build();
    }
}
