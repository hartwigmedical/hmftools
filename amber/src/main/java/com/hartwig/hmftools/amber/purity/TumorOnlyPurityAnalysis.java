package com.hartwig.hmftools.amber.purity;

import static java.lang.String.format;

import static com.google.common.collect.Range.open;
import static com.hartwig.hmftools.amber.AmberConfig.AMB_LOGGER;
import static com.hartwig.hmftools.amber.AmberConstants.CHR_ARM_AUC_UPPER_BOUND_FOR_COPY_NUMBER_EVENTS;
import static com.hartwig.hmftools.amber.AmberConstants.CHR_ARM_LOWER_BOUND_FOR_CONTAMINATION;
import static com.hartwig.hmftools.amber.AmberConstants.HET_VAF_LOWER_BOUND;
import static com.hartwig.hmftools.amber.AmberConstants.HET_VAF_UPPER_BOUND;
import static com.hartwig.hmftools.amber.AmberConstants.HOMOZYGOUS_PROPORTION_LOWER_BOUND_FOR_CONTAMINATION;
import static com.hartwig.hmftools.amber.AmberConstants.HOMOZYGOUS_PROPORTION_UPPER_BOUND_FOR_CONTAMINATION;
import static com.hartwig.hmftools.amber.AmberConstants.MIN_CUTOFF;
import static com.hartwig.hmftools.amber.AmberConstants.MUTATION_AUC_LOWER_BOUND_FOR_CONTAMINATION;

import java.util.ArrayList;
import java.util.List;
import java.util.function.Function;

import com.google.common.collect.ListMultimap;
import com.google.common.collect.Range;
import com.hartwig.hmftools.amber.PositionEvidence;
import com.hartwig.hmftools.common.amber.AmberSite;
import com.hartwig.hmftools.common.genome.chromosome.Chromosome;
import com.hartwig.hmftools.common.genome.refgenome.RefGenomeVersion;
import com.hartwig.hmftools.common.immune.ImmuneRegions;
import com.hartwig.hmftools.common.region.ChrBaseRegion;
import com.hartwig.hmftools.common.segmentation.ChrArm;
import com.hartwig.hmftools.common.segmentation.ChrArmLocator;
import com.hartwig.hmftools.common.utils.Doubles;

public class TumorOnlyPurityAnalysis
{
    public static final Range<Double> HOM_PROP_RANGE =
            open(HOMOZYGOUS_PROPORTION_LOWER_BOUND_FOR_CONTAMINATION, HOMOZYGOUS_PROPORTION_UPPER_BOUND_FOR_CONTAMINATION);
    private final List<CandidatePeak> mCopyNumberPeaks = new ArrayList<>();
    private final List<CandidatePeak> mContaminationPeaks = new ArrayList<>();

    static Function<PositionEvidence, CanonicalSnvType> snvTypeClassifier()
    {
        return evidence ->
        {
            if(evidence.vaf() <= HET_VAF_UPPER_BOUND)
            {
                return CanonicalSnvType.type(evidence.Ref, evidence.Alt);
            }
            else
            {
                return CanonicalSnvType.type(evidence.Alt, evidence.Ref);
            }
        };
    }

    static Function<PositionEvidence, ChrArm> chrArmClassifier(final RefGenomeVersion genomeVersion)
    {
        ChrArmLocator locator = ChrArmLocator.defaultLocator(genomeVersion);
        return evidence -> locator.map(evidence.chromosome(), evidence.position());
    }

    public TumorOnlyPurityAnalysis(
            final List<PositionEvidence> evidence,
            final ListMultimap<Chromosome, AmberSite> amberSites,
            final PurityAnalysisConfig config)
    {
        AMB_LOGGER.debug(format("evaluating noise floor with %d evidence points", evidence.size()));
        List<PositionEvidence> filteredEvidence = filterOutExcludedRegions(evidence, config);
        AMB_LOGGER.debug(format("immune region data removed leaving %d evidence points", filteredEvidence.size()));
        List<PositionEvidence> hetVariants = getBaselineHetVariants(filteredEvidence);
        GnomadFrequencySupplier frequencySupplier = new DefaultGnomadFrequencySupplier(amberSites, config.refGenomeVersion());
        GnomadFrequencyAnalysis baselineGnomadAnalysis = new GnomadFrequencyAnalysis(frequencySupplier);
        double baselineHetGnomadFrequency = baselineGnomadAnalysis.getMeanFrequency(hetVariants);
        AMB_LOGGER.debug(format("baseline het gnomad frequency: %.3f", baselineHetGnomadFrequency));

        AucCalculator<PositionEvidence, ChrArm> chrArmAucCalculator =
                new AucCalculator<>(chrArmClassifier(config.refGenomeVersion()), hetVariants);
        AucCalculator<PositionEvidence, CanonicalSnvType> snvTypeAucCalculator = new AucCalculator<>(snvTypeClassifier(), hetVariants);
        for(CandidatePeakEvaluationResult evaluationResult : new PeakSearch(filteredEvidence, config.threads()).peaks())
        {
            CandidatePeak peak = evaluationResult.candidatePeak();
            PeakGnomadFrequenciesChecker gnomadChecker = new PeakGnomadFrequenciesChecker(peak);
            if(!gnomadChecker.checkGnomadFrequencies(frequencySupplier, baselineHetGnomadFrequency))
            {
                AMB_LOGGER.debug(format("peak %.3f failed Gnomad frequency check", peak.vaf()));
                continue;
            }
            double homozygousProportion = peak.homozygousProportion();
            double chrArmAuc = chrArmAucCalculator.calculateAuc(peak.allCapturedPoints());
            double mutationAuc = snvTypeAucCalculator.calculateAuc(peak.allCapturedPoints());
            AMB_LOGGER.debug(format("peak %.3f, homProportion: %.3f, chrArmAUC: %.3f, mutationAUC: %.3f", peak.vaf(), homozygousProportion, chrArmAuc, mutationAuc));
            if(HOM_PROP_RANGE.contains(homozygousProportion))
            {
                if(mutationAuc > MUTATION_AUC_LOWER_BOUND_FOR_CONTAMINATION
                        && chrArmAuc > CHR_ARM_LOWER_BOUND_FOR_CONTAMINATION)
                {
                    AMB_LOGGER.debug(format("peak %.3f is a contamination peak", peak.vaf()));
                    mContaminationPeaks.add(peak);
                    continue;
                }
            }
            if(chrArmAuc < CHR_ARM_AUC_UPPER_BOUND_FOR_COPY_NUMBER_EVENTS
                    || mutationAuc > MUTATION_AUC_LOWER_BOUND_FOR_CONTAMINATION && chrArmAuc < CHR_ARM_LOWER_BOUND_FOR_CONTAMINATION)
            {
                AMB_LOGGER.debug(format("peak %.3f is a copy number peak", peak.vaf()));
                mCopyNumberPeaks.add(peak);
            }
        }
    }

    private static List<PositionEvidence> filterOutExcludedRegions(final List<PositionEvidence> evidence, final PurityAnalysisConfig config)
    {
        List<ChrBaseRegion> immuneRegions = getExcludedImmuneRegions(config.refGenomeVersion());
        RegionsFilter immuneRegionsFilter = new RegionsFilter(immuneRegions);
        return immuneRegionsFilter.filter(evidence);
    }

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
        double minCopyNumberPeak = mCopyNumberPeaks.stream().map(CandidatePeak::vaf).min(Double::compare).orElse(0.0);
        if(Doubles.isZero(minCopyNumberPeak))
        {
            return MIN_CUTOFF;
        }
        return Math.min(MIN_CUTOFF, minCopyNumberPeak / 3.0);
    }

    public List<CandidatePeak> contaminationPeaks()
    {
        return mContaminationPeaks;
    }

    public List<CandidatePeak> copyNumberPeaks()
    {
        return mCopyNumberPeaks;
    }

    private static List<PositionEvidence> getBaselineHetVariants(final List<PositionEvidence> evidence)
    {
        List<PositionEvidence> baselineHetPositions = new ArrayList<>();
        for(PositionEvidence pe : evidence)
        {
            final double vaf = pe.symmetricVaf();
            if(vaf > HET_VAF_LOWER_BOUND && vaf < HET_VAF_UPPER_BOUND)
            {
                baselineHetPositions.add(pe);
            }
        }
        return baselineHetPositions;
    }
}
