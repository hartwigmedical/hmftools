package com.hartwig.hmftools.amber.purity;

import static java.lang.String.format;

import static com.hartwig.hmftools.amber.AmberConfig.AMB_LOGGER;

import java.util.ArrayList;
import java.util.List;

import com.google.common.collect.ListMultimap;
import com.hartwig.hmftools.amber.PositionEvidence;
import com.hartwig.hmftools.common.amber.AmberSite;
import com.hartwig.hmftools.common.genome.chromosome.Chromosome;
import com.hartwig.hmftools.common.genome.position.GenomePosition;
import com.hartwig.hmftools.common.genome.refgenome.RefGenomeVersion;
import com.hartwig.hmftools.common.segmentation.ChrArmLocator;

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

    public TumorOnlyNoiseFloorAnalysis(List<PositionEvidence> evidence,
            RefGenomeVersion refGenomeVersion,
            ListMultimap<Chromosome, AmberSite> amberSites)
    {
        PeakSearch search = new PeakSearch(evidence);
        GnomadFrequencySupplier frequencySupplier = new DefaultGnomadFrequencySupplier(amberSites, refGenomeVersion);
        double baselineHetGnomadFrequency = getBaselineHetGnomadFrequency(evidence, frequencySupplier);
        AMB_LOGGER.debug(format("Baseline Het Gnomad Frequency: %.3f", baselineHetGnomadFrequency));

        for(VafLevelEvaluationResult evaluationResult : search.peaks())
        {
            VafLevel peak = evaluationResult.Vaf();
            PeakGnomadFrequenciesChecker gnomadChecker = new PeakGnomadFrequenciesChecker(peak);
            boolean ok = gnomadChecker.checkGnomadFrequencies(frequencySupplier, baselineHetGnomadFrequency);
            if(!ok)
            {
                AMB_LOGGER.debug(format("Peak %.3f failed Gnomad frequency check", peak.vaf()));
                continue;
            }
            double homozygousProportion = peak.homozygousProportion();
            double chrArmAuc = peak.perArmConsistencyFactor(ChrArmLocator.defaultLocator(refGenomeVersion));
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
