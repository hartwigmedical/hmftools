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
import com.hartwig.hmftools.common.region.BaseRegion;
import com.hartwig.hmftools.common.region.ChrBaseRegion;
import com.hartwig.hmftools.common.utils.pcf.PCFPosition;
import com.hartwig.hmftools.purple.AmberData;
import com.hartwig.hmftools.purple.CobaltData;

import org.apache.commons.math3.stat.descriptive.moment.StandardDeviation;

public class ChimerismDetection
{
    private final AmberData mAmberData;
    private final CobaltData mCobaltData;
    private final RefGenomeVersion mRefGenomeVersion;

    private final List<AmberPcfRegion> mAmberPcfRegions;

    // results
    private boolean mIsDetected;

    public ChimerismDetection(final AmberData amberData, final CobaltData cobaltData, final RefGenomeVersion refGenomeVersion)
    {
        mAmberData = amberData;
        mCobaltData = cobaltData;
        mRefGenomeVersion = refGenomeVersion;

        mAmberPcfRegions = Lists.newArrayList();
        mIsDetected = false;
    }

    public void run()
    {
        detectChrimerism();
    }

    public boolean isDetected() { return mIsDetected; }

    private void detectChrimerism()
    {
        for(Chromosome chromosome : mAmberData.TumorSegments.keySet())
        {
            String chrStr = mRefGenomeVersion.versionedChromosome(chromosome.toString());
            List<BaseRegion> pcfRegions = mAmberData.PcfRegions.get(chrStr);

            List<AmberBAF> amberBAFs = mAmberData.ChromosomeBafs.get(chromosome).stream().collect(Collectors.toList());

            int amberBafIndex = 0;

            for(BaseRegion pcfRegion : pcfRegions)
            {
                if(pcfRegion.baseLength() <= 1)
                    continue;

                AmberPcfRegion amberPcfRegion = new AmberPcfRegion(chrStr, pcfRegion.start(), pcfRegion.end());

                while(amberBafIndex < amberBAFs.size())
                {
                    AmberBAF amberBAF = amberBAFs.get(amberBafIndex);

                    if(amberBAF.Position > pcfRegion.end())
                        break;

                    if(amberBAF.Position >= pcfRegion.start())
                    {
                        amberPcfRegion.AmberBAFs.add(amberBAF);
                    }

                    ++amberBafIndex;
                }

                if(amberPcfRegion.bafCount() >= 2)
                {
                    mAmberPcfRegions.add(amberPcfRegion);
                }
            }
        }

        // calculate a weighted average from all regions
        double weightedBafTotal = 0;
        double countsTotal = 0;

        for(AmberPcfRegion pcfRegion : mAmberPcfRegions)
        {
            weightedBafTotal += pcfRegion.bafStandardDeviation() * pcfRegion.bafCount();
            countsTotal += pcfRegion.bafCount();
        }

        double sampleAverage = weightedBafTotal / countsTotal;

        mIsDetected = sampleAverage > CHIMERISM_SAMPLE_CUTOFF;
    }

    private class AmberPcfRegion extends ChrBaseRegion
    {
        public final List<AmberBAF> AmberBAFs;

        private Double mStandardDeviation;

        public AmberPcfRegion(final String chromosome, final int posStart, final int posEnd)
        {
            super(chromosome, posStart, posEnd);
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
