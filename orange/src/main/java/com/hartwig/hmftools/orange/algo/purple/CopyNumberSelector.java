package com.hartwig.hmftools.orange.algo.purple;

import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.stream.Collectors;

import com.google.common.collect.Maps;
import com.google.common.collect.Sets;
import com.hartwig.hmftools.common.driver.panel.DriverGene;
import com.hartwig.hmftools.common.genome.chromosome.HumanChromosome;
import com.hartwig.hmftools.common.purple.GeneCopyNumber;
import com.hartwig.hmftools.datamodel.purple.CopyNumberInterpretation;
import com.hartwig.hmftools.datamodel.purple.ImmutablePurpleGainDeletion;
import com.hartwig.hmftools.datamodel.purple.PurpleGainDeletion;

import static com.hartwig.hmftools.orange.OrangeApplication.LOGGER;

import com.google.common.collect.Lists;
import org.jetbrains.annotations.NotNull;

final class CopyNumberSelector
{
    @NotNull
    public static List<PurpleGainDeletion> selectNearReportableSomaticGains(@NotNull List<GeneCopyNumber> allGeneCopyNumbers, double ploidy,
            @NotNull List<PurpleGainDeletion> reportableGainsLosses, @NotNull List<DriverGene> driverGenes)
    {
        List<PurpleGainDeletion> nearReportableSomaticGains = Lists.newArrayList();
        Set<String> ampDriverGenes = selectAmpDriverGenes(driverGenes);
        for(GeneCopyNumber geneCopyNumber : allGeneCopyNumbers)
        {
            if(ampDriverGenes.contains(geneCopyNumber.geneName()))
            {
                double relativeMinCopyNumber = geneCopyNumber.minCopyNumber() / ploidy;
                double relativeMaxCopyNumber = geneCopyNumber.maxCopyNumber() / ploidy;
                if(relativeMinCopyNumber > 2.5 && relativeMaxCopyNumber < 3)
                {
                    nearReportableSomaticGains.add(toFullGain(geneCopyNumber));
                }
            }
        }

        // Check in case official amp have changed.
        Set<String> reportableGenes = Sets.newHashSet();
        for(PurpleGainDeletion reportable : reportableGainsLosses)
        {
            reportableGenes.add(reportable.gene());
        }

        for(PurpleGainDeletion gain : nearReportableSomaticGains)
        {
            if(reportableGenes.contains(gain.gene()))
            {
                LOGGER.warn("Gene {} is selected to be near-reportable but has already been reported!", gain.gene());
            }
        }

        return nearReportableSomaticGains;
    }

    @NotNull
    public static List<PurpleGainDeletion> selectInterestingUnreportedGainsLosses(@NotNull List<PurpleGainDeletion> allGainsLosses,
            @NotNull List<PurpleGainDeletion> reportableGainsLosses)
    {
        List<PurpleGainDeletion> unreportedGainDels = selectUnreportedGainsLosses(allGainsLosses, reportableGainsLosses);

        List<PurpleGainDeletion> interestingUnreportedGainsLosses = Lists.newArrayList();
        interestingUnreportedGainsLosses.addAll(selectInterestingGains(unreportedGainDels));
        interestingUnreportedGainsLosses.addAll(selectInterestingLosses(unreportedGainDels, reportableGainsLosses));
        return interestingUnreportedGainsLosses;
    }

    @NotNull
    private static Set<String> selectAmpDriverGenes(@NotNull List<DriverGene> driverGenes)
    {
        Set<String> ampGenes = Sets.newHashSet();
        for(DriverGene driverGene : driverGenes)
        {
            if(driverGene.reportAmplification())
            {
                ampGenes.add(driverGene.gene());
            }
        }
        return ampGenes;
    }

    @NotNull
    private static PurpleGainDeletion toFullGain(@NotNull GeneCopyNumber geneCopyNumber)
    {
        return ImmutablePurpleGainDeletion.builder()
                .chromosome(geneCopyNumber.chromosome())
                .chromosomeBand(geneCopyNumber.chromosomeBand())
                .gene(geneCopyNumber.geneName())
                .transcript(geneCopyNumber.transName())
                .isCanonical(geneCopyNumber.isCanonical())
                .interpretation(CopyNumberInterpretation.FULL_GAIN)
                .minCopies(Math.max(0, geneCopyNumber.minCopyNumber()))
                .maxCopies(Math.max(0, geneCopyNumber.maxCopyNumber()))
                .build();
    }

    @NotNull
    private static List<PurpleGainDeletion> selectUnreportedGainsLosses(@NotNull List<PurpleGainDeletion> allGainsLosses,
            @NotNull List<PurpleGainDeletion> reportableGainsLosses)
    {
        List<PurpleGainDeletion> unreportedGainsLosses = Lists.newArrayList();
        for(PurpleGainDeletion gainDel : allGainsLosses)
        {
            if(!reportableGainsLosses.contains(gainDel))
            {
                unreportedGainsLosses.add(gainDel);
            }
        }
        return unreportedGainsLosses;
    }

    @NotNull
    private static List<PurpleGainDeletion> selectInterestingGains(@NotNull List<PurpleGainDeletion> unreportedgainDels)
    {
        List<PurpleGainDeletion> unreportedFullGains = unreportedgainDels.stream()
                .filter(gainDel -> gainDel.interpretation() == CopyNumberInterpretation.FULL_GAIN)
                .collect(Collectors.toList());

        Map<CopyNumberKey, PurpleGainDeletion> bestGainPerLocation = Maps.newHashMap();
        for(PurpleGainDeletion gain : unreportedFullGains)
        {
            CopyNumberKey key = new CopyNumberKey(gain.chromosome(), gain.chromosomeBand());
            PurpleGainDeletion bestGain = bestGainPerLocation.get(key);
            if(bestGain == null)
            {
                bestGainPerLocation.put(key, gain);
            }
            else
            {
                if(gain.minCopies() > bestGain.minCopies())
                {
                    bestGainPerLocation.put(key, gain);
                }
            }
        }

        return Lists.newArrayList(bestGainPerLocation.values().iterator());
    }

    @NotNull
    private static List<PurpleGainDeletion> selectInterestingLosses(@NotNull List<PurpleGainDeletion> unreportedGainsLosses,
            @NotNull List<PurpleGainDeletion> reportableGainsLosses)
    {
        List<PurpleGainDeletion> unreportedLosses = unreportedGainsLosses.stream()
                .filter(gainDel -> gainDel.interpretation() == CopyNumberInterpretation.PARTIAL_DEL
                        || gainDel.interpretation() == CopyNumberInterpretation.FULL_DEL)
                .collect(Collectors.toList());

        List<PurpleGainDeletion> reportableLosses = reportableGainsLosses.stream()
                .filter(gainDel -> gainDel.interpretation() == CopyNumberInterpretation.PARTIAL_DEL
                        || gainDel.interpretation() == CopyNumberInterpretation.FULL_DEL)
                .collect(Collectors.toList());

        List<PurpleGainDeletion> lossesAutosomes = Lists.newArrayList();
        for(PurpleGainDeletion loss : unreportedLosses)
        {
            if(HumanChromosome.fromString(loss.chromosome()).isAutosome())
            {
                if(!locusPresent(reportableLosses, loss.chromosome(), loss.chromosomeBand()))
                {
                    lossesAutosomes.add(loss);
                }
            }
        }

        Map<CopyNumberKey, PurpleGainDeletion> bestLossPerLocation = Maps.newHashMap();
        for(PurpleGainDeletion loss : lossesAutosomes)
        {
            CopyNumberKey key = new CopyNumberKey(loss.chromosome(), loss.chromosomeBand());
            PurpleGainDeletion bestLoss = bestLossPerLocation.get(key);
            if(bestLoss == null)
            {
                bestLossPerLocation.put(key, loss);
            }
            else
            {
                boolean pickOtherWhenEqual = bestLoss.gene().compareTo(loss.gene()) <= 0;
                if(pickOtherWhenEqual)
                {
                    bestLossPerLocation.put(key, loss);
                }
            }
        }

        return Lists.newArrayList(bestLossPerLocation.values().iterator());
    }

    private static boolean locusPresent(@NotNull List<PurpleGainDeletion> gainsLosses, @NotNull String chromosome,
            @NotNull String chromosomeBand)
    {
        for(PurpleGainDeletion gainDel : gainsLosses)
        {
            if(gainDel.chromosome().equals(chromosome) && gainDel.chromosomeBand().equals(chromosomeBand))
            {
                return true;
            }
        }

        return false;
    }
}
