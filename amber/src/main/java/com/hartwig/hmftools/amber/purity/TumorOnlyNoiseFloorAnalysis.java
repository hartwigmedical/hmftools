package com.hartwig.hmftools.amber.purity;

import static java.lang.String.format;

import static com.hartwig.hmftools.amber.AmberConfig.AMB_LOGGER;

import java.util.ArrayList;
import java.util.List;

import com.google.common.collect.ListMultimap;
import com.hartwig.hmftools.amber.AmberConfig;
import com.hartwig.hmftools.amber.PositionEvidence;
import com.hartwig.hmftools.amber.PositionEvidenceFile;
import com.hartwig.hmftools.common.amber.AmberSite;
import com.hartwig.hmftools.common.genome.chromosome.Chromosome;
import com.hartwig.hmftools.common.genome.position.GenomePosition;
import com.hartwig.hmftools.common.genome.refgenome.RefGenomeVersion;
import com.hartwig.hmftools.common.immune.ImmuneRegions;
import com.hartwig.hmftools.common.region.ChrBaseRegion;
import com.hartwig.hmftools.common.segmentation.ChrArmLocator;

import org.jetbrains.annotations.NotNull;

public class TumorOnlyNoiseFloorAnalysis
{
    private static final double HET_VAF_LOWER_BOUND = 0.35;
    private static final double HET_VAF_UPPER_BOUND = 0.65;
    private static final double HOMOZYGOUS_PROPORTION_LOWER_BOUND_FOR_CONTAMINATION = 0.25;
    private static final double HOMOZYGOUS_PROPORTION_UPPER_BOUND_FOR_CONTAMINATION = 0.75;
    private static final double MUTATION_AUC_LOWER_BOUND_FOR_CONTAMINATION = 0.825;
    private static final double CHR_ARM_LOWER_BOUND_FOR_CONTAMINATION = 0.3;
    private static final double CHR_ARM_AUC_UPPER_BOUND_FOR_COPY_NUMBER_EVENTS = 0.1;
    private static final double MIN_CUTOFF = 0.04;
    private List<VafLevel> CopyNumberPeaks = new ArrayList<>();
    private List<VafLevel> ContaminationPeaks = new ArrayList<>();

    public TumorOnlyNoiseFloorAnalysis(
            final List<PositionEvidence> evidence,
            final ListMultimap<Chromosome, AmberSite> amberSites,
            final AmberConfig config) throws Exception
    {
        String rawDataFileName = PositionEvidenceFile.generateTumorDataFilename(config.OutputDir, config.TumorId);
        PositionEvidenceFile.write(rawDataFileName, evidence);
        AMB_LOGGER.debug(format("Evaluating noise floor with %d evidence points", evidence.size()));
        List<ChrBaseRegion> immuneRegions = getExcludedImmuneRegions(config.RefGenVersion);
        RegionsFilter immuneRegionsFilter = new RegionsFilter(immuneRegions);
        List<PositionEvidence> filteredEvidence = immuneRegionsFilter.filter(evidence);
        AMB_LOGGER.debug(format("Immune region data removed leaving %d evidence points", filteredEvidence.size()));
        PeakSearch search = new PeakSearch(filteredEvidence, config.Threads);
        GnomadFrequencySupplier frequencySupplier = new DefaultGnomadFrequencySupplier(amberSites, config.RefGenVersion);
        double baselineHetGnomadFrequency = getBaselineHetGnomadFrequency(filteredEvidence, frequencySupplier);
        AMB_LOGGER.debug(format("Baseline Het Gnomad Frequency: %.3f", baselineHetGnomadFrequency));

        for(VafLevelEvaluationResult evaluationResult : search.peaks())
        {
            VafLevel peak = evaluationResult.Vaf();
            PeakGnomadFrequenciesChecker gnomadChecker = new PeakGnomadFrequenciesChecker(peak);
            boolean okByGnomadFrequencies = gnomadChecker.checkGnomadFrequencies(frequencySupplier, baselineHetGnomadFrequency);
            if(!okByGnomadFrequencies)
            {
                AMB_LOGGER.debug(format("Peak %.3f failed Gnomad frequency check", peak.vaf()));
                continue;
            }
            double homozygousProportion = peak.homozygousProportion();
            double chrArmAuc = peak.perArmConsistencyFactor(ChrArmLocator.defaultLocator(config.RefGenVersion));
            double mutationAuc = peak.perMutationTypeConsistencyFactor();
            AMB_LOGGER.debug(format("Peak %.3f, homProportion: %.3f, chrArmAUC: %.3f, mutationAUC: %.3f", peak.vaf(), homozygousProportion, chrArmAuc, mutationAuc));
            if(homozygousProportion > HOMOZYGOUS_PROPORTION_LOWER_BOUND_FOR_CONTAMINATION
                    && homozygousProportion < HOMOZYGOUS_PROPORTION_UPPER_BOUND_FOR_CONTAMINATION)
            {
                if(mutationAuc > MUTATION_AUC_LOWER_BOUND_FOR_CONTAMINATION
                        && chrArmAuc > CHR_ARM_LOWER_BOUND_FOR_CONTAMINATION)
                {
                    AMB_LOGGER.debug(format("Peak %.3f is a contamination peak", peak.vaf()));
                    ContaminationPeaks.add(peak);
                    continue;
                }
            }
            if(chrArmAuc < CHR_ARM_AUC_UPPER_BOUND_FOR_COPY_NUMBER_EVENTS
                    || mutationAuc > MUTATION_AUC_LOWER_BOUND_FOR_CONTAMINATION && chrArmAuc < CHR_ARM_LOWER_BOUND_FOR_CONTAMINATION)
            {
                AMB_LOGGER.debug(format("Peak %.3f is a copy number peak", peak.vaf()));
                CopyNumberPeaks.add(peak);
            }
        }
    }

    @NotNull
    private static List<ChrBaseRegion> getExcludedImmuneRegions(final RefGenomeVersion refGenomeVersion)
    {
        List<ChrBaseRegion> immuneRegions = new ArrayList<>();
        immuneRegions.addAll(ImmuneRegions.getHlaRegions(refGenomeVersion));
        immuneRegions.addAll(ImmuneRegions.getIgRegions(refGenomeVersion));
        immuneRegions.addAll(ImmuneRegions.getTrRegions(refGenomeVersion));
        return immuneRegions;
    }

    public double cutoff()
    {
        double minCopyNumberPeak = CopyNumberPeaks.stream().map(VafLevel::vaf).min(Double::compare).orElse(0.0);
        if(minCopyNumberPeak == 0.0)
        {
            return MIN_CUTOFF;
        }
        return minCopyNumberPeak / 3.0;
    }

    public List<VafLevel> copyNumberPeaks()
    {
        return CopyNumberPeaks;
    }

    public List<VafLevel> contaminationPeaks()
    {
        return ContaminationPeaks;
    }

    private static double getBaselineHetGnomadFrequency(final List<PositionEvidence> evidence,
            final GnomadFrequencySupplier frequencySupplier)
    {
        GnomadFrequencyAnalysis baselineGnomadAnalysis = new GnomadFrequencyAnalysis(frequencySupplier);
        List<GenomePosition> baselineHetPositions = new ArrayList<>();
        for(PositionEvidence pe : evidence)
        {
            final double vaf = pe.symmetricVaf();
            if(vaf > HET_VAF_LOWER_BOUND && vaf < HET_VAF_UPPER_BOUND)
            {
                baselineHetPositions.add(pe);
            }
        }
        return baselineGnomadAnalysis.getMeanFrequency(baselineHetPositions);
    }
}
