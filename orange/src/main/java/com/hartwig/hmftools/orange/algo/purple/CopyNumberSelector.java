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
import com.hartwig.hmftools.datamodel.purple.ImmutablePurpleGainLoss;
import com.hartwig.hmftools.datamodel.purple.PurpleGainDel;

import static com.hartwig.hmftools.orange.OrangeApplication.LOGGER;

import com.google.common.collect.Lists;
import org.jetbrains.annotations.NotNull;

final class CopyNumberSelector
{
    @NotNull
    public static List<PurpleGainDel> selectNearReportableSomaticGains(@NotNull List<GeneCopyNumber> allGeneCopyNumbers, double ploidy,
            @NotNull List<PurpleGainDel> reportableGainsLosses, @NotNull List<DriverGene> driverGenes)
    {
        List<PurpleGainDel> nearReportableSomaticGains = Lists.newArrayList();
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
        for(PurpleGainDel reportable : reportableGainsLosses)
        {
            reportableGenes.add(reportable.gene());
        }

        for(PurpleGainDel gain : nearReportableSomaticGains)
        {
            if(reportableGenes.contains(gain.gene()))
            {
                LOGGER.warn("Gene {} is selected to be near-reportable but has already been reported!", gain.gene());
            }
        }

        return nearReportableSomaticGains;
    }

    @NotNull
    public static List<PurpleGainDel> selectInterestingUnreportedGainsLosses(@NotNull List<PurpleGainDel> allGainsLosses,
            @NotNull List<PurpleGainDel> reportableGainsLosses)
    {
        List<PurpleGainDel> unreportedGainLosses = selectUnreportedGainsLosses(allGainsLosses, reportableGainsLosses);

        List<PurpleGainDel> interestingUnreportedGainsLosses = Lists.newArrayList();
        interestingUnreportedGainsLosses.addAll(selectInterestingGains(unreportedGainLosses));
        interestingUnreportedGainsLosses.addAll(selectInterestingLosses(unreportedGainLosses, reportableGainsLosses));
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
    private static PurpleGainDel toFullGain(@NotNull GeneCopyNumber geneCopyNumber)
    {
        return ImmutablePurpleGainLoss.builder()
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
    private static List<PurpleGainDel> selectUnreportedGainsLosses(@NotNull List<PurpleGainDel> allGainsLosses,
            @NotNull List<PurpleGainDel> reportableGainsLosses)
    {
        List<PurpleGainDel> unreportedGainsLosses = Lists.newArrayList();
        for(PurpleGainDel gainLoss : allGainsLosses)
        {
            if(!reportableGainsLosses.contains(gainLoss))
            {
                unreportedGainsLosses.add(gainLoss);
            }
        }
        return unreportedGainsLosses;
    }

    @NotNull
    private static List<PurpleGainDel> selectInterestingGains(@NotNull List<PurpleGainDel> unreportedGainLosses)
    {
        List<PurpleGainDel> unreportedFullGains = unreportedGainLosses.stream()
                .filter(gainLoss -> gainLoss.interpretation() == CopyNumberInterpretation.FULL_GAIN)
                .collect(Collectors.toList());

        Map<CopyNumberKey, PurpleGainDel> bestGainPerLocation = Maps.newHashMap();
        for(PurpleGainDel gain : unreportedFullGains)
        {
            CopyNumberKey key = new CopyNumberKey(gain.chromosome(), gain.chromosomeBand());
            PurpleGainDel bestGain = bestGainPerLocation.get(key);
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
    private static List<PurpleGainDel> selectInterestingLosses(@NotNull List<PurpleGainDel> unreportedGainsLosses,
            @NotNull List<PurpleGainDel> reportableGainsLosses)
    {
        List<PurpleGainDel> unreportedLosses = unreportedGainsLosses.stream()
                .filter(gainLoss -> gainLoss.interpretation() == CopyNumberInterpretation.PARTIAL_DEL
                        || gainLoss.interpretation() == CopyNumberInterpretation.FULL_DEL)
                .collect(Collectors.toList());

        List<PurpleGainDel> reportableLosses = reportableGainsLosses.stream()
                .filter(gainLoss -> gainLoss.interpretation() == CopyNumberInterpretation.PARTIAL_DEL
                        || gainLoss.interpretation() == CopyNumberInterpretation.FULL_DEL)
                .collect(Collectors.toList());

        List<PurpleGainDel> lossesAutosomes = Lists.newArrayList();
        for(PurpleGainDel loss : unreportedLosses)
        {
            if(HumanChromosome.fromString(loss.chromosome()).isAutosome())
            {
                if(!locusPresent(reportableLosses, loss.chromosome(), loss.chromosomeBand()))
                {
                    lossesAutosomes.add(loss);
                }
            }
        }

        Map<CopyNumberKey, PurpleGainDel> bestLossPerLocation = Maps.newHashMap();
        for(PurpleGainDel loss : lossesAutosomes)
        {
            CopyNumberKey key = new CopyNumberKey(loss.chromosome(), loss.chromosomeBand());
            PurpleGainDel bestLoss = bestLossPerLocation.get(key);
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

    private static boolean locusPresent(@NotNull List<PurpleGainDel> gainsLosses, @NotNull String chromosome,
            @NotNull String chromosomeBand)
    {
        for(PurpleGainDel gainLoss : gainsLosses)
        {
            if(gainLoss.chromosome().equals(chromosome) && gainLoss.chromosomeBand().equals(chromosomeBand))
            {
                return true;
            }
        }

        return false;
    }
}
