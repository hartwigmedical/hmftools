package com.hartwig.hmftools.redux.bqr;

import static java.lang.Math.log10;
import static java.lang.Math.max;
import static java.lang.Math.min;

import static com.hartwig.hmftools.common.codon.Nucleotides.DNA_BASE_BYTES;
import static com.hartwig.hmftools.common.redux.BaseQualAdjustment.probabilityToPhredQual;
import static com.hartwig.hmftools.redux.ReduxConfig.RD_LOGGER;
import static com.hartwig.hmftools.redux.ReduxConstants.BQR_CHR_END_BUFFER;
import static com.hartwig.hmftools.redux.ReduxConstants.BQR_SAMPLE_SIZE;

import java.util.Collection;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.stream.Collectors;

import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import com.google.common.collect.Sets;
import com.hartwig.hmftools.common.genome.chromosome.HumanChromosome;
import com.hartwig.hmftools.common.genome.refgenome.RefGenomeCoordinates;
import com.hartwig.hmftools.common.redux.BqrFile;
import com.hartwig.hmftools.common.redux.BqrKey;
import com.hartwig.hmftools.common.redux.BqrRecord;
import com.hartwig.hmftools.common.region.ChrBaseRegion;
import com.hartwig.hmftools.common.utils.RExecutor;
import com.hartwig.hmftools.redux.ReduxConfig;

public class BaseQualRecalibration
{
    private final ReduxConfig mConfig;

    private final boolean mEnabled;
    private final BaseQualityResults mResults;

    private final List<ChrBaseRegion> mRegions;

    public BaseQualRecalibration(final ReduxConfig config)
    {
        mConfig = config;
        mResults = new BaseQualityResults();

        mRegions = Lists.newArrayList();

        if(mConfig.BQR.Enabled)
        {
            mEnabled = true;
            formRegions();
        }
        else
        {
            mEnabled = false;
        }
    }

    public boolean enabled() { return mEnabled; }
    public BaseQualityResults results() { return mResults; }
    public List<ChrBaseRegion> regions() { return mRegions; }

    private void formRegions()
    {
        // mimic the regions created in Sage, or take them all if running over full BAM
        if(!mConfig.SpecificChrRegions.Regions.isEmpty())
        {
            for(ChrBaseRegion region : mConfig.SpecificChrRegions.Regions)
            {
                mRegions.add(new ChrBaseRegion(region.Chromosome, region.start(), region.end()));
            }

            return;
        }

        // form regions from 2MB per chromosome
        RefGenomeCoordinates refGenomeCoordinates = RefGenomeCoordinates.refGenomeCoordinates(mConfig.RefGenVersion);

        for(HumanChromosome chromosome : HumanChromosome.values())
        {
            if(chromosome.isAllosome())
                continue;

            String chrStr = mConfig.RefGenVersion.versionedChromosome(chromosome.toString());

            if(mConfig.SpecificChrRegions.excludeChromosome(chrStr))
                continue;

            List<ChrBaseRegion> chrRegionTasks = Lists.newArrayList();

            int chrLength = refGenomeCoordinates.length(chrStr);

            int chrRegionStart, chrRegionEnd;

            if(mConfig.BQR.UseAllRegions)
            {
                chrRegionStart = BQR_CHR_END_BUFFER;
                chrRegionEnd = chrLength - BQR_CHR_END_BUFFER;
            }
            else
            {
                chrRegionEnd = chrLength - BQR_CHR_END_BUFFER;
                chrRegionStart = chrRegionEnd - BQR_SAMPLE_SIZE;
            }

            /*
            int start = chrRegionStart;

            while(start < chrRegionEnd)
            {
                chrRegionTasks.add(new ChrBaseRegion(chrStr, start, start + BQR_REGION_SIZE - 1));
                start += BQR_REGION_SIZE;
            }
            */

            chrRegionTasks.add(new ChrBaseRegion(chrStr, chrRegionStart, chrRegionEnd));

            mRegions.addAll(chrRegionTasks);
        }
    }

    public void finalise()
    {
        RD_LOGGER.debug("writing base quality recalibration, readsUsed({}) altsFiltered({})",
                mResults.totalReadsUsed(), mResults.totalAltsFiltered());

        String sampleId = mConfig.SampleId;;

        // merge results for this sample across all regions
        Map<BqrKey,Integer> allQualityCounts = mResults.getCombinedQualityCounts();

        List<BqrRecord> records = convertToRecords(allQualityCounts);

        // write results to file
        writeSampleData(sampleId, records);
    }

    private static final byte NO_BASE_OR_QUAL = 1; // value is irrelevant, just to complete a map entry

    public static List<BqrRecord> convertToRecords(final Map<BqrKey,Integer> allQualityCounts)
    {
        List<BqrRecord> result = Lists.newArrayList();

        // collect the ref counts
        Map<BqrKey,Integer> refCountMap = allQualityCounts.entrySet().stream()
                .filter(x -> x.getKey().Ref == x.getKey().Alt)
                .collect(Collectors.toMap(x -> x.getKey(), x -> x.getValue()));

        // make a map of (per-type) trinuc totals across all entries
        Map<BqrKey,Integer> triNucMap = Maps.newHashMap();

        for(Map.Entry<BqrKey,Integer> entry : allQualityCounts.entrySet())
        {
            BqrKey key = entry.getKey();
            int count = entry.getValue();

            if(key.Quality == 0)
                continue;

            BqrKey triNucKey = new BqrKey(NO_BASE_OR_QUAL, NO_BASE_OR_QUAL, key.TrinucleotideContext, NO_BASE_OR_QUAL, key.ReadType);

            Integer triNucCount = triNucMap.get(triNucKey);
            triNucMap.put(triNucKey, triNucCount != null ? triNucCount + count : count);
        }

        Set<BqrKey> syntheticAltKeys = Sets.newHashSet();

        for(Map.Entry<BqrKey,Integer> entry : allQualityCounts.entrySet())
        {
            BqrKey key = entry.getKey();

            if(key.Quality == 0)
                continue;

            double recalibratedQual = 0;

            BqrKey refKey = new BqrKey(key.Ref, key.Ref, key.TrinucleotideContext, key.Quality, key.ReadType);

            if(key.Alt == key.Ref)
            {
                recalibratedQual = key.Quality;
            }
            else
            {
                recalibratedQual = calcRecalibratedQual(key, entry.getValue(), refCountMap, triNucMap);
            }

            result.add(new BqrRecord(key, entry.getValue(), recalibratedQual));

            // add alt entries for any which sample has no results
            for(int i = 0; i < DNA_BASE_BYTES.length; ++i)
            {
                byte alt = DNA_BASE_BYTES[i];

                if(alt == refKey.Ref)
                    continue;

                BqrKey altKey = new BqrKey(key.Ref, alt, key.TrinucleotideContext, key.Quality, key.ReadType);

                if(!allQualityCounts.containsKey(altKey) && !syntheticAltKeys.contains(altKey))
                {
                    syntheticAltKeys.add(altKey);

                    double recalQualMin = calcRecalibratedQual(altKey, 1, refCountMap, triNucMap);
                    double syntheticQual = max(recalQualMin + 10 * log10(2), key.Quality);

                    result.add(new BqrRecord(altKey, 0, syntheticQual));
                }
            }
        }

        return result;
    }

    private static double calcRecalibratedQual(
            final BqrKey key, final int observedCount, final Map<BqrKey,Integer> refCountMap, final Map<BqrKey,Integer> triNucMap)
    {
        // example: to calc qual of T>C errors at a GTA trinucleotide, with qual 30
        // error_rate = (# all GCA sites / # all GTA sites) * (# qual 30 T>Cs GTA sites) / (# qual 30 T>Cs GTA sites + # qual 30 C refs GCA sites)

        // or in general terms:
        // errorRate = (# alt TN sites / # ref TN sites) * (# alts at ref TN sites) / (# alts at ref TN sites + # refs at alt TN sites)

        byte[] altTriNucContext = new byte[] { key.TrinucleotideContext[0], key.Alt, key.TrinucleotideContext[2] };

        BqrKey altKey = new BqrKey(key.Alt, key.Alt, altTriNucContext, key.Quality, key.ReadType);
        int altRefCount = refCountMap.getOrDefault(altKey, 0);

        BqrKey altTriNucKey = new BqrKey(NO_BASE_OR_QUAL, NO_BASE_OR_QUAL, altTriNucContext, NO_BASE_OR_QUAL, key.ReadType);
        int altTriNucCount = triNucMap.getOrDefault(altTriNucKey, 0);

        BqrKey refTriNucKey = new BqrKey(NO_BASE_OR_QUAL, NO_BASE_OR_QUAL, key.TrinucleotideContext, NO_BASE_OR_QUAL, key.ReadType);
        int refTriNucCount = triNucMap.getOrDefault(refTriNucKey, 0);

        double triNucRate = refTriNucCount > 0 ? altTriNucCount / (double)refTriNucCount : 0;

        if(triNucRate == 0 || (observedCount + altRefCount) == 0)
            return 0;

        double calcProbability = min(triNucRate * observedCount / (observedCount + altRefCount), 1.0);
        return probabilityToPhredQual(calcProbability);
    }

    private void writeSampleData(final String sampleId, final Collection<BqrRecord> records)
    {
        try
        {
            String tsvFile = BqrFile.generateFilename(mConfig.OutputDir, sampleId);

            BqrFile.write(tsvFile, records.stream().collect(Collectors.toList()));

            if(mConfig.BQR.WritePlot)
            {
                RExecutor.executeFromClasspath("r/baseQualityRecalibrationPlot.R", tsvFile);
            }
        }
        catch(Exception e)
        {
            RD_LOGGER.error(" sample({}) failed to write base recalibration: {}", sampleId, e.toString());
        }
    }
}
