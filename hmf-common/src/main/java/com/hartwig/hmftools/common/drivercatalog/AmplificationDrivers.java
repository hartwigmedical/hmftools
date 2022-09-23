package com.hartwig.hmftools.common.drivercatalog;

import static com.hartwig.hmftools.common.drivercatalog.DriverCatalogFactory.createCopyNumberDriver;

import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.function.Predicate;
import java.util.stream.Collectors;

import com.google.common.collect.Lists;
import com.google.common.collect.Sets;
import com.hartwig.hmftools.common.drivercatalog.panel.DriverGene;
import com.hartwig.hmftools.common.drivercatalog.panel.DriverGenePanel;
import com.hartwig.hmftools.common.genome.chromosome.HumanChromosome;
import com.hartwig.hmftools.common.purple.PurpleQCStatus;
import com.hartwig.hmftools.common.purple.GeneCopyNumber;
import com.hartwig.hmftools.common.purple.SegmentSupport;
import com.hartwig.hmftools.common.utils.Doubles;

public class AmplificationDrivers
{
    private static final double MIN_COPY_NUMBER_RELATIVE_INCREASE = 3;

    public static final double MAX_COPY_NUMBER_DEL = 0.5;

    private static final List<String> TARGET_REGIONS_PARTIAL_AMP_GENES = Lists.newArrayList(
            "BRAF", "EGFR", "CTNNB1", "CBL", "MET", "ALK", "PDGFRA");

    private final Set<PurpleQCStatus> mQcStatus;
    private final Map<String, DriverGene> mAmplificationTargets;

    public AmplificationDrivers(final Set<PurpleQCStatus> qcStatus, final DriverGenePanel panel)
    {
        mQcStatus = qcStatus;
        mAmplificationTargets = panel.amplificationTargets().stream().collect(Collectors.toMap(DriverGene::gene, x -> x));
    }

    public List<DriverCatalog> amplifications(final double ploidy, final List<GeneCopyNumber> geneCopyNumbers, boolean isTargetRegions)
    {
        List<DriverCatalog> result = Lists.newArrayList();

        boolean isHighCopyNoise = mQcStatus.contains(PurpleQCStatus.WARN_HIGH_COPY_NUMBER_NOISE);

        for(GeneCopyNumber geneCopyNumber : geneCopyNumbers)
        {
            if(!mAmplificationTargets.containsKey(geneCopyNumber.geneName()))
                continue;

            if(isHighCopyNoise && !supportedByOneSV(geneCopyNumber))
                continue;

            if(geneCopyNumber.minCopyNumber() / ploidy > MIN_COPY_NUMBER_RELATIVE_INCREASE)
            {
                result.add(createAmpDriver(geneCopyNumber));
                continue;
            }

            if(!isTargetRegions || TARGET_REGIONS_PARTIAL_AMP_GENES.contains(geneCopyNumber.geneName()))
            {
                if(geneCopyNumber.maxCopyNumber() / ploidy > MIN_COPY_NUMBER_RELATIVE_INCREASE)
                {
                    result.add(createPartialAmpDriver(geneCopyNumber));
                }
            }
        }

        return result;
    }

    static boolean supportedByOneSV(final GeneCopyNumber geneCopyNumber)
    {
        return geneCopyNumber.minRegionStartSupport().isSV() || geneCopyNumber.minRegionEndSupport().isSV();
    }

    private DriverCatalog createAmpDriver(final GeneCopyNumber geneCopyNumber)
    {
        DriverGene driverGene = mAmplificationTargets.get(geneCopyNumber.geneName());
        return createCopyNumberDriver(driverGene.likelihoodType(), DriverType.AMP, LikelihoodMethod.AMP, false, geneCopyNumber);
    }

    private DriverCatalog createPartialAmpDriver(final GeneCopyNumber geneCopyNumber)
    {
        DriverGene driverGene = mAmplificationTargets.get(geneCopyNumber.geneName());

        return createCopyNumberDriver(driverGene.likelihoodType(), DriverType.PARTIAL_AMP, LikelihoodMethod.AMP, false, geneCopyNumber);
    }
}
