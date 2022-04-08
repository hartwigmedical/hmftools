package com.hartwig.hmftools.cobalt.targeted

import com.google.common.collect.ArrayListMultimap
import com.google.common.collect.ListMultimap
import com.google.common.collect.Sets
import com.hartwig.hmftools.cobalt.Chromosome
import com.hartwig.hmftools.cobalt.CobaltConstants
import com.hartwig.hmftools.cobalt.ratio.RatioBuilder
import com.hartwig.hmftools.common.cobalt.ImmutableReadRatio
import com.hartwig.hmftools.common.cobalt.ReadRatio
import com.hartwig.hmftools.common.genome.position.GenomePosition
import com.hartwig.hmftools.common.genome.window.Window
import com.hartwig.hmftools.common.utils.Doubles
import org.apache.logging.log4j.Level
import org.apache.logging.log4j.LogManager
import kotlin.math.roundToInt

class TargetedRatioBuilder(
    targetRegions: List<GenomePosition>,
    targetEnrichmentRatios: Map<GenomePosition, Double>,
    rawRatios: ListMultimap<Chromosome, ReadRatio>) : RatioBuilder
{
    private val mOnTargetRatios = ArrayListMultimap.create<Chromosome, ReadRatio>()
    private val mOffTargetRatios = ArrayListMultimap.create<Chromosome, ReadRatio>()
    private val mCombinedRatios = ArrayListMultimap.create<Chromosome, ReadRatio>()

    // we use on target ratios only for now
    override fun ratios(): ArrayListMultimap<Chromosome, ReadRatio>
    {
        return mOnTargetRatios
    }

    fun onTargetRatios(): ArrayListMultimap<Chromosome, ReadRatio>
    {
        return mOnTargetRatios
    }

    fun offTargetRatios(): ArrayListMultimap<Chromosome, ReadRatio>
    {
        return mOffTargetRatios
    }

    init
    {
        // calculate median for sample
        populateOnTargetRatios(targetEnrichmentRatios, rawRatios)
        populateOffTargetRatios(rawRatios, CobaltConstants.OFF_TARGET_WINDOW_SIZE, targetRegions)
        sLogger.info("{} on target GC ratios, {} off target GC ratios", mOnTargetRatios.size(), mOffTargetRatios.size())
        populateCombinedRatios(mOnTargetRatios, mOffTargetRatios)
    }

    private fun populateCombinedRatios(ratios1: ArrayListMultimap<Chromosome, ReadRatio>,
                                       ratios2: ArrayListMultimap<Chromosome, ReadRatio>)
    {
        mCombinedRatios.clear()
        val chromosomes = Sets.newIdentityHashSet<Chromosome>()
        chromosomes.addAll(ratios1.keySet())
        chromosomes.addAll(ratios2.keySet())
        for (chromosome in chromosomes)
        {
            var i1 = 0
            var i2 = 0
            val list1 = ratios1[chromosome]
            val list2 = ratios2[chromosome]
            while (i1 < list1.size && i2 < list2.size)
            {
                val readRatio1 = list1[i1]
                val readRatio2 = list2[i2]
                if (readRatio1.position() < readRatio2.position())
                {
                    mCombinedRatios.put(chromosome, readRatio1)
                    ++i1
                }
                else if (readRatio1.position() > readRatio2.position())
                {
                    mCombinedRatios.put(chromosome, readRatio2)
                    ++i2
                }
                else
                {
                    // shouldn't be here
                    sLogger.error(
                        "on target and off target ratio at same position: {}:{}",
                        readRatio1.chromosome(), readRatio1.position())
                    throw RuntimeException("on target and off target ratio at same position")
                }
            }

            // put the rest in
            while (i1 < list1.size)
            {
                mCombinedRatios.put(chromosome, list1[i1])
                ++i1
            }
            while (i2 < list2.size)
            {
                mCombinedRatios.put(chromosome, list2[i2])
                ++i2
            }
        }
    }

    private fun populateOnTargetRatios(
        targetRelativeEnrichment: Map<GenomePosition, Double>,
        rawRatios: ListMultimap<Chromosome, ReadRatio>)
    {
        // find all the ratios that are inside the target enriched regions
        // we filter out all the regions with 0 gc normalised ratios, as they do not actually
        // correctly reflect the amount of enrichment, and also very rare
        val targetRegionsGcRatios: List<Double> = rawRatios.entries()
            .filter { (_, readRatio): Map.Entry<Chromosome, ReadRatio> ->
                readRatio.ratio() > 0.0 && targetRelativeEnrichment.containsKey(readRatio)
            }
            .map { (_, value): Map.Entry<Chromosome, ReadRatio> -> value.ratio() }

        var targetRegionGcRatioMedian = 1.0

        if (targetRegionsGcRatios.isNotEmpty())
        {
            targetRegionGcRatioMedian = Doubles.median(targetRegionsGcRatios)
        }
        else
        {
            sLogger.warn("target region gc ratios is empty")
        }

        sLogger.printf(Level.INFO, "targeted mode GC ratio median: %.3f", targetRegionGcRatioMedian)
        mOnTargetRatios.clear()
        for ((key, value) in rawRatios.entries())
        {
            val relativeEnrichment = targetRelativeEnrichment[value] ?: continue
            val enrichmentAdjRatio = value.ratio() / targetRegionGcRatioMedian / relativeEnrichment
            sLogger.debug(
                "{}:{} relative enrichment: {}, on target ratio: {}",
                value.chromosome(),
                value.position(),
                relativeEnrichment,
                enrichmentAdjRatio)
            mOnTargetRatios.put(key, ImmutableReadRatio.builder().from(value).ratio(enrichmentAdjRatio).build())
        }
    }

    // we create a pan window ratio by taking the median count of super windows that combine multiple windows
    private fun populateOffTargetRatios(
        rawRatios: ListMultimap<Chromosome, ReadRatio>, offTargetWindowSize: Int,
        targetRegions: List<GenomePosition>)
    {
        val unnormalizedRatios = ArrayListMultimap.create<Chromosome, ReadRatio>()
        val window = Window(offTargetWindowSize)
        for (chromosome in rawRatios.keySet())
        {
            var currentWindowStart: Int = -1
            val windowGcRatios = ArrayList<Double>()

            // we need this to make sure we get consistent chromosome name (1 vs chr1)
            val chromosomeStr: String = rawRatios[chromosome][0].chromosome()
            for (readRatio in rawRatios[chromosome])
            {
                // todo: make sure this is sorted
                val windowStart: Int = window.start(readRatio.position())
                if (windowStart != currentWindowStart)
                {
                    if (currentWindowStart != -1)
                    {
                        val unnormalizedRatio = unnormalizedOffTargetRatio(
                            offTargetWindowSize,
                            chromosomeStr,
                            currentWindowStart,
                            windowGcRatios,
                            targetRegions)
                        if (unnormalizedRatio != null)
                        {
                            unnormalizedRatios.put(chromosome, unnormalizedRatio)
                        }
                    }
                    currentWindowStart = windowStart
                    windowGcRatios.clear()
                }
                if (readRatio.ratio() >= 0.0)
                {
                    windowGcRatios.add(readRatio.ratio())
                }
            }
            val unnormalizedRatio = unnormalizedOffTargetRatio(
                offTargetWindowSize,
                chromosomeStr,
                currentWindowStart,
                windowGcRatios,
                targetRegions)
            if (unnormalizedRatio != null)
            {
                unnormalizedRatios.put(chromosome, unnormalizedRatio)
            }
        }

        // now we want to normalise all of those off target gc ratios by the median
        val median: Double = Doubles.median(unnormalizedRatios.values().map({ obj: ReadRatio -> obj.ratio() }))
        sLogger.debug("normalizing {} off target windows ratio by median: {}", unnormalizedRatios.size(), median)
        mOffTargetRatios.clear()
        for ((key, value) in unnormalizedRatios.entries())
        {
            val normalizedRatio = value.ratio() / median
            mOffTargetRatios.put(
                key,
                ImmutableReadRatio.builder().from(value).ratio(normalizedRatio).build())
        }
    }

    companion object
    {
        private val sLogger = LogManager.getLogger(TargetedRatioBuilder::class.java)
        private fun unnormalizedOffTargetRatio(offTargetWindowSize: Int,
                                               chromosome: String,
                                               windowStart: Int,
                                               windowGcRatios: List<Double>,
                                               targetRegions: List<GenomePosition>): ReadRatio?
        {
            val minNumGcRatios =
                (offTargetWindowSize.toDouble() / CobaltConstants.WINDOW_SIZE * CobaltConstants.MIN_OFF_TARGET_WINDOW_RATIO).roundToInt()
            if (windowGcRatios.size < minNumGcRatios)
            {
                // if we don't have enough sub windows with valid values then we skip this
                return null
            }
            val windowEnd = windowStart + offTargetWindowSize - 1

            // check that this window does not contain any target regions
            for (targetRegion in targetRegions)
            {
                if (windowStart <= targetRegion.position() && windowEnd > targetRegion.position())
                {
                    // this window contains a target region
                    return null
                }
            }
            val median = Doubles.median(windowGcRatios)

            // the window position is the middle
            val windowMid = windowStart + offTargetWindowSize / 2
            sLogger.debug(
                "off target window: {}:{} ({} - {}), num sub windows: {}, median: {}",
                chromosome, windowMid,
                windowStart, windowEnd, windowGcRatios.size, median)
            if (!median.isNaN())
            {
                return ImmutableReadRatio.builder().chromosome(chromosome).position(windowMid).ratio(median).build()
            }
            return null
        }
    }
}