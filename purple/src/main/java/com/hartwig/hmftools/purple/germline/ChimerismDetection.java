package com.hartwig.hmftools.purple.germline;

import static java.lang.String.format;

import static com.hartwig.hmftools.purple.PurpleConstants.CHIMERISM_KDE_BANDWIDTH;
import static com.hartwig.hmftools.purple.PurpleConstants.CHIMERISM_MIN_BAF_COUNT;
import static com.hartwig.hmftools.purple.PurpleConstants.CHIMERISM_SAMPLE_CUTOFF;
import static com.hartwig.hmftools.purple.PurpleUtils.PPL_LOGGER;

import java.util.List;
import java.util.Map;
import java.util.stream.Collectors;
import java.util.stream.DoubleStream;
import java.util.stream.IntStream;

import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import com.hartwig.hmftools.common.amber.AmberBAF;
import com.hartwig.hmftools.common.genome.chromosome.Chromosome;
import com.hartwig.hmftools.common.genome.refgenome.RefGenomeVersion;
import com.hartwig.hmftools.common.purple.GermlineStatus;
import com.hartwig.hmftools.common.region.ChrBaseRegion;
import com.hartwig.hmftools.common.utils.Doubles;
import com.hartwig.hmftools.common.utils.kde.KernelEstimator;
import com.hartwig.hmftools.purple.AmberData;
import com.hartwig.hmftools.purple.CobaltData;
import com.hartwig.hmftools.purple.region.ObservedRegion;

import org.apache.commons.math3.stat.descriptive.moment.StandardDeviation;

public class ChimerismDetection
{
    private final AmberData mAmberData;
    private final CobaltData mCobaltData;
    private List<ObservedRegion> mFilteredRegions;
    private final RefGenomeVersion mRefGenomeVersion;

    private final List<RegionBafData> mRegionBafData;

    // results
    private double mChimerismLevel;

    public ChimerismDetection(
            final AmberData amberData, final CobaltData cobaltData, List<ObservedRegion> observedRegions, final RefGenomeVersion refGenomeVersion)
    {
        mAmberData = amberData;
        mCobaltData = cobaltData;
        mRefGenomeVersion = refGenomeVersion;

        mRegionBafData = Lists.newArrayList();

        mFilteredRegions = observedRegions.stream()
                .filter(x -> x.germlineStatus() == GermlineStatus.DIPLOID)
                .filter(x -> x.bafCount() >= CHIMERISM_MIN_BAF_COUNT)
                .collect(Collectors.toList());

        mChimerismLevel = 0;
    }

    public boolean isDetected() { return mChimerismLevel > CHIMERISM_SAMPLE_CUTOFF; }

    public void run()
    {
        detectChrimerism();

        applyFit();
    }

    private void detectChrimerism()
    {
        Map<String,List<ObservedRegion>> chrFilteredRegions = Maps.newHashMap();

        for(ObservedRegion region : mFilteredRegions)
        {
            List<ObservedRegion> regions = chrFilteredRegions.get(region.chromosome());

            if(regions == null)
            {
                regions = Lists.newArrayList();
                chrFilteredRegions.put(region.chromosome(), regions);
            }

            regions.add(region);
        }

        for(Chromosome chromosome : mAmberData.TumorSegments.keySet())
        {
            String chrStr = mRefGenomeVersion.versionedChromosome(chromosome.toString());

            List<ObservedRegion> filteredRegions = chrFilteredRegions.get(chrStr);

            if(filteredRegions == null)
                continue;

            List<AmberBAF> amberBAFs = mAmberData.ChromosomeBafs.get(chromosome).stream().collect(Collectors.toList());

            int amberBafIndex = 0;

            for(ObservedRegion region : filteredRegions)
            {
                RegionBafData regionBafData = new RegionBafData(region);

                while(amberBafIndex < amberBAFs.size())
                {
                    AmberBAF amberBAF = amberBAFs.get(amberBafIndex);

                    if(amberBAF.Position > regionBafData.end())
                        break;

                    if(amberBAF.Position >= regionBafData.start())
                    {
                        regionBafData.AmberBAFs.add(amberBAF);
                    }

                    ++amberBafIndex;
                }

                if(regionBafData.bafCount() >= 2)
                {
                    mRegionBafData.add(regionBafData);
                }
            }
        }

        // calculate a weighted average from all regions
        double weightedBafTotal = 0;
        double countsTotal = 0;

        for(RegionBafData pcfRegion : mRegionBafData)
        {
            weightedBafTotal += pcfRegion.bafStandardDeviation() * pcfRegion.bafCount();
            countsTotal += pcfRegion.bafCount();
        }

        mChimerismLevel = weightedBafTotal / countsTotal;
    }

    private void applyFit()
    {
        final KernelEstimator estimator = new KernelEstimator(0.001, CHIMERISM_KDE_BANDWIDTH);

        List<Double> bafValues = Lists.newArrayList();

        for(ObservedRegion region : mFilteredRegions)
        {
            // TODO: get cobalt data, calculate a BAF
            double baf = 0;

            bafValues.add(baf);
            estimator.addValue(baf, 1.0);
        }

        double[] bafs = IntStream.rangeClosed(0, 51).mapToDouble(x -> x / 100d).toArray();
        double[] densities = DoubleStream.of(bafs).map(estimator::getProbability).toArray();

        for(int i = 1; i < densities.length - 1; i++)
        {
            double density = densities[i];
            if(Doubles.greaterThan(density, densities[i - 1]) && Doubles.greaterThan(density, densities[i + 1]))
            {
                double baf = bafs[i];

                int peakCount = count(baf, bafValues);

                PPL_LOGGER.debug(format("discovered peak: count(%d) baf(%.3f)", peakCount, baf));
            }
        }
    }

    private static int count(double peak, final List<Double> values)
    {
        // TODO: is a buffer around the peak required?
        return (int) values.stream().filter(vaf -> between(peak, vaf - 0.015, vaf + 0.015)).count();
    }

    private static boolean between(double victim, double min, double max)
    {
        return Doubles.greaterOrEqual(victim, min) && Doubles.lessOrEqual(victim, max);
    }

    private class RegionBafData extends ChrBaseRegion
    {
        public final ObservedRegion Region;
        public final List<AmberBAF> AmberBAFs;

        private Double mStandardDeviation;

        public RegionBafData(final ObservedRegion region)
        {
            super(region.chromosome(), region.start(), region.end());
            Region = region;

            AmberBAFs = Lists.newArrayList();
            mStandardDeviation = null;
        }

        public double bafStandardDeviation()
        {
            if(mStandardDeviation != null)
                return mStandardDeviation;

            double[] doubles = new double[AmberBAFs.size()];

            for(int i = 0; i < AmberBAFs.size(); ++i)
            {
                doubles[i] = AmberBAFs.get(i).tumorModifiedBAF();
            }

            mStandardDeviation = new StandardDeviation().evaluate(doubles);
            return mStandardDeviation;
        }

        public int bafCount() { return AmberBAFs.size(); }

        public String toString() { return format("region(%s) points(%d)", super.toString(), AmberBAFs.size()); }
    }
}
