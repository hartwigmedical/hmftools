package com.hartwig.hmftools.purple.germline;

import static java.lang.String.format;

import static com.hartwig.hmftools.purple.PurpleConstants.CHIMERISM_SAMPLE_CUTOFF;

import java.util.List;
import java.util.Map;
import java.util.stream.Collectors;

import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import com.hartwig.hmftools.common.amber.AmberBAF;
import com.hartwig.hmftools.common.genome.chromosome.Chromosome;
import com.hartwig.hmftools.common.genome.refgenome.RefGenomeVersion;
import com.hartwig.hmftools.common.purple.GermlineStatus;
import com.hartwig.hmftools.common.region.ChrBaseRegion;
import com.hartwig.hmftools.purple.AmberData;
import com.hartwig.hmftools.purple.region.ObservedRegion;

import org.apache.commons.math3.stat.descriptive.moment.StandardDeviation;

public class ChimerismDetection
{
    private final AmberData mAmberData;
    private List<ObservedRegion> mObservedRegions;
    private final RefGenomeVersion mRefGenomeVersion;

    private final List<RegionBafData> mRegionBafData;

    // results
    private boolean mIsDetected;

    public ChimerismDetection(final AmberData amberData, List<ObservedRegion> observedRegions, final RefGenomeVersion refGenomeVersion)
    {
        mAmberData = amberData;
        mObservedRegions = observedRegions;
        mRefGenomeVersion = refGenomeVersion;

        mRegionBafData = Lists.newArrayList();
        mIsDetected = false;
    }

    public void run()
    {
        detectChrimerism();
    }

    public boolean isDetected() { return mIsDetected; }

    private void detectChrimerism()
    {
        Map<String,List<ObservedRegion>> chrFilteredRegions = Maps.newHashMap();

        for(ObservedRegion region : mObservedRegions)
        {
            if(region.germlineStatus() != GermlineStatus.DIPLOID)
                continue;

            if(region.bafCount() < 2)
                continue;

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

        double sampleAverage = weightedBafTotal / countsTotal;

        mIsDetected = sampleAverage > CHIMERISM_SAMPLE_CUTOFF;
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
