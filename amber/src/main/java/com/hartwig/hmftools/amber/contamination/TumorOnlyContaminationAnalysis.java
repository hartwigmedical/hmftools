package com.hartwig.hmftools.amber.contamination;

import static java.lang.String.format;

import static com.hartwig.hmftools.amber.AmberConfig.AMB_LOGGER;
import static com.hartwig.hmftools.amber.contamination.PerClassVafConsistencyChecker.calculateConfirmationFactor;
import static com.hartwig.hmftools.common.segmentation.ChrArmLocator.defaultLocator;

import java.util.ArrayList;
import java.util.Collection;
import java.util.List;
import java.util.SortedSet;
import java.util.TreeSet;

import com.hartwig.hmftools.amber.PositionEvidence;
import com.hartwig.hmftools.amber.VafReading;
import com.hartwig.hmftools.amber.TumorBAF;
import com.hartwig.hmftools.common.amber.BaseDepthData;
import com.hartwig.hmftools.common.genome.refgenome.RefGenomeVersion;
import com.hartwig.hmftools.common.segmentation.ChrArm;

public class TumorOnlyContaminationAnalysis
{
    interface VariantFrequencySupplier
    {
        double minorAlleleFrequencyInPopulation(TumorBAF tumorBAF);
    }

    List<VafReading> VafReadings = new ArrayList<>();
    List<TumorContamination> ContaminationPoints = new ArrayList<>();
    private final double ContaminationEstimate;
    private final VafConsistencyCheckResult<ChrArm> Check;

    public TumorOnlyContaminationAnalysis(final Collection<TumorBAF> tumorBafs, final RefGenomeVersion refGenVersion,
            VariantFrequencySupplier frequencySupplier)
    {
        SortedSet<VafReading> readDepthData = new TreeSet<>();
        SortedSet<TumorContamination> contaminationSortedSet = new TreeSet<>();
        for(TumorBAF tumorBAF : tumorBafs)
        {
            PositionEvidence x = tumorBAF.TumorEvidence;
            if(useEvidenceForContaminationDetection(x))
            {
                double minorAlleleFrequency = frequencySupplier.minorAlleleFrequencyInPopulation(tumorBAF);
                readDepthData.add(new VafReading(x.Chromosome, x.Position, x.ReadDepth, x.RefSupport, x.AltSupport, minorAlleleFrequency));
                BaseDepthData tumorBdd = new BaseDepthData(x.Ref, x.Alt, x.ReadDepth, x.IndelCount, x.RefSupport, x.AltSupport);
                contaminationSortedSet.add(new TumorContamination(x.Chromosome, x.Position, null, tumorBdd));
            }
        }
        VafReadings.addAll(readDepthData);
        ContaminationPoints.addAll(contaminationSortedSet);
        if(ContaminationPoints.size() < 0.01 * tumorBafs.size())
        {
            AMB_LOGGER.debug("Too few contamination points to estimate contamination.");
            ContaminationEstimate = 0.0;
            Check = null;
        }
        else
        {
            SearchGrid searchGrid = new SearchGrid();
            SearchGrid.Calculator calculator = new CumulativeProbabilityContaminationCost(VafReadings);
            SearchGrid.ValueScore score = searchGrid.findBestValue(calculator);
            if(score == null)
            {
                ContaminationEstimate = 0.0;
                Check = null;
            }
            else
            {
                ContaminationEstimate = score.value();
                AMB_LOGGER.debug(format("Initial tumor contamination estimate: %.3f", ContaminationEstimate));
                Check = calculateConfirmationFactor(defaultLocator(refGenVersion), ContaminationEstimate, VafReadings);
                AMB_LOGGER.debug(format("Contamination check: %.3f,vaf weight: %d, total weight: %d", Check.unevenDistributionCost(), Check.totalWeightInBand(), Check.totalWeightAcrossAllVafValues()));
            }
        }
    }

    private static boolean useEvidenceForContaminationDetection(final PositionEvidence x)
    {
        if(x.ReadDepth <= 7)
        {
            return false;
        }
        double v = x.vaf();
        return v > 0.004 && v < 0.4 || v > 0.6 && v < 0.996;
    }

    public List<TumorContamination> getContaminationPoints()
    {
        return ContaminationPoints;
    }

    public double getApproximateContamination()
    {
        if(Check != null && Check.unevenDistributionCost() > 0.25)
        {
            return 0.0;
        }
        return ContaminationEstimate;
    }

    public List<CategoryEvidence<ChrArm>> getArmEvidence()
    {
        return Check.categoryEvidence();
    }
}
