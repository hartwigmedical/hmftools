package com.hartwig.hmftools.common.driver;

import static com.hartwig.hmftools.common.driver.DriverCatalogFactory.createCopyNumberDriver;

import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.stream.Collectors;

import com.google.common.collect.Lists;
import com.google.common.collect.Sets;
import com.hartwig.hmftools.common.driver.panel.DriverGene;
import com.hartwig.hmftools.common.driver.panel.DriverGenePanel;
import com.hartwig.hmftools.common.genome.chromosome.HumanChromosome;
import com.hartwig.hmftools.common.purple.GeneCopyNumber;
import com.hartwig.hmftools.common.purple.PurpleQCStatus;
import com.hartwig.hmftools.common.purple.SegmentSupport;
import com.hartwig.hmftools.common.utils.Doubles;

public class DeletionDrivers
{
    public static final double MAX_COPY_NUMBER_DEL = 0.5;

    private static final int SHORT_DEL_LENGTH = 10_000_000;
    private static final double SINGLE_DEPTH_GC_MIN = 0.35;
    private static final double SINGLE_DEPTH_GC_MAX = 0.6;

    private static final Set<SegmentSupport> MERE = Sets.newHashSet(SegmentSupport.CENTROMERE, SegmentSupport.TELOMERE);

    private final Set<PurpleQCStatus> mQcStatus;
    private final Map<String, DriverGene> mDeletionTargets;

    public DeletionDrivers(final Set<PurpleQCStatus> qcStatus, final DriverGenePanel panel)
    {
        mQcStatus = qcStatus;
        mDeletionTargets = panel.deletionTargets().stream().collect(Collectors.toMap(DriverGene::gene, x -> x));
    }

    public static int deletedGenes(final List<GeneCopyNumber> geneCopyNumbers)
    {
        return (int) geneCopyNumbers.stream()
                .filter(x -> !HumanChromosome.fromString(x.chromosome()).equals(HumanChromosome._Y)
                        && Doubles.lessThan(x.minCopyNumber(), MAX_COPY_NUMBER_DEL))
                .count();
    }

    public List<DriverCatalog> deletions(final List<GeneCopyNumber> geneCopyNumbers, boolean isTargetRegions)
    {
        List<DriverCatalog> drivers = Lists.newArrayList();

        boolean checkQcStatus = mQcStatus.contains(PurpleQCStatus.WARN_DELETED_GENES)
                || mQcStatus.contains(PurpleQCStatus.WARN_HIGH_COPY_NUMBER_NOISE);

        for(GeneCopyNumber geneCopyNumber : geneCopyNumbers)
        {
            if(geneCopyNumber.minCopyNumber() >= MAX_COPY_NUMBER_DEL)
                continue;

            if(!mDeletionTargets.containsKey(geneCopyNumber.geneName()))
                continue;

            boolean hasFullSvSupport = supportedByTwoSVs(geneCopyNumber);
            boolean isShort = geneCopyNumber.minRegionBases() < SHORT_DEL_LENGTH;

            if(isTargetRegions)
            {
                if(!hasFullSvSupport)
                {
                    if(checkQcStatus)
                        continue;

                    if(geneCopyNumber.depthWindowCount() == 1)
                    {
                        if(geneCopyNumber.gcContent() < SINGLE_DEPTH_GC_MIN || geneCopyNumber.gcContent() > SINGLE_DEPTH_GC_MAX)
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

            drivers.add(createDelDriver(geneCopyNumber));
        }

        return drivers;
    }

    static boolean supportedByTwoSVs(final GeneCopyNumber geneCopyNumber)
    {
        return geneCopyNumber.minRegionStartSupport().isSV() && geneCopyNumber.minRegionEndSupport().isSV();
    }

    static boolean supportedByOneSVAndMere(final GeneCopyNumber geneCopyNumber)
    {
        if(MERE.contains(geneCopyNumber.minRegionStartSupport()) && geneCopyNumber.minRegionEndSupport().isSV())
        {
            return true;
        }

        if(MERE.contains(geneCopyNumber.minRegionEndSupport()) && geneCopyNumber.minRegionStartSupport().isSV())
        {
            return true;
        }

        return false;
    }

    private DriverCatalog createDelDriver(final GeneCopyNumber geneCopyNumber)
    {
        DriverGene driverGene = mDeletionTargets.get(geneCopyNumber.geneName());
        return createCopyNumberDriver(driverGene.likelihoodType(), DriverType.DEL, LikelihoodMethod.DEL, true, geneCopyNumber);
    }
}
