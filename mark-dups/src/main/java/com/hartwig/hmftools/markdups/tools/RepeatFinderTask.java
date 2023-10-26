package com.hartwig.hmftools.markdups.tools;

import static com.hartwig.hmftools.common.genome.refgenome.RefGenomeSource.loadRefGenome;
import static com.hartwig.hmftools.common.genome.refgenome.RefGenomeVersion.V37;
import static com.hartwig.hmftools.markdups.tools.RepeatFinderConfig.MD_LOGGER;

import java.io.BufferedWriter;
import java.util.List;
import java.util.concurrent.Callable;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.genome.refgenome.RefGenomeCoordinates;
import com.hartwig.hmftools.common.genome.refgenome.RefGenomeSource;
import com.hartwig.hmftools.common.region.BaseRegion;
import com.hartwig.hmftools.common.region.ChrBaseRegion;
import com.hartwig.hmftools.common.utils.PerformanceCounter;

import org.jetbrains.annotations.Nullable;

public class RepeatFinderTask implements Callable
{
    private final RepeatFinderConfig mConfig;
    private final String mChromosome;
    private final BufferedWriter mWriter;
    private boolean mBedAnnotations;
    private final List<BaseRegion> mChrBedRegions;

    private final RefGenomeSource mRefGenome;
    private final RefGenomeCoordinates mRefGenomeCoords;

    private final PerformanceCounter mPerfCounter;
    private int mSingleNuclRepeatCount;
    private int mDoubleNuclRepeatCount;
    private int mTotalSingleNuclRepeatBases;
    private int mTotalDoubleNuclRepeatBases;

    public RepeatFinderTask(final RepeatFinderConfig config, final String chromosome, final BufferedWriter writer)
    {
        this(config, chromosome, writer, null);
        mBedAnnotations = false;
    }

    public RepeatFinderTask(final RepeatFinderConfig config, final String chromosome, final BufferedWriter writer,
            @Nullable final List<BaseRegion> chrBedRegions)
    {
        mConfig = config;
        mChromosome = chromosome;
        mWriter = writer;
        mBedAnnotations = true;
        mChrBedRegions = (chrBedRegions == null) ? Lists.newArrayList() : chrBedRegions;

        mRefGenome = loadRefGenome(config.RefGenome);
        mRefGenomeCoords = mConfig.RefGenVersion == V37 ? RefGenomeCoordinates.COORDS_37 : RefGenomeCoordinates.COORDS_38;

        mPerfCounter = new PerformanceCounter(String.format("RepeatFinderTask Chr(%s)", mChromosome));
        mSingleNuclRepeatCount = 0;
        mDoubleNuclRepeatCount = 0;
        mTotalSingleNuclRepeatBases = 0;
        mTotalDoubleNuclRepeatBases = 0;
    }

    @Override
    public Long call()
    {
        MD_LOGGER.info("RepeatFinderTask started for chromosome {}", mChromosome);

        mPerfCounter.start();
        processChromosome();
        mPerfCounter.stop();

        mPerfCounter.logStats();
        MD_LOGGER.info(
                "RepeatFinderTask for chromosome {} finished, {} single nucleotide repeats found that cover a total of {} bases, {} double nucleotide repeats found that cover a total of {} bases",
                mChromosome, mSingleNuclRepeatCount, mTotalSingleNuclRepeatBases, mDoubleNuclRepeatCount, mTotalDoubleNuclRepeatBases);

        return (long) 0;
    }

    /**
     * <p>Finds single and double nucleotide base repeats whose length exceed `RepeatFinderConfig.MinRepeatBases`.</p>
     * <p>
     * Note that we treat a single nucleotide repeat as a double nucleotide repeat. Note that for a double nucleotide repeat to make
     * sense, it should be at least four bases long. Whereas, a two base repeat of a single nucleotide could be considered a repeat.
     * However, treating a single nucleotide repeat as a double nucleotide repeat means that this will only find single nucleotide repeats
     * of length at least four. However, we are interested in much longer repeats than this, so this isn't an issue.
     * </p>
     */
    private void processChromosome()
    {
        final int chromosomeLength = mRefGenomeCoords.length(mChromosome);
        final String refBases = mRefGenome.getBaseString(mChromosome, 1, chromosomeLength);

        final List<Character> history = Lists.newArrayList();
        int currentRepeatLength = 0;
        int chrBedRegionIdx = 0;
        boolean shiftChrBedRegionIdx = false;
        for(int i = 0; i < chromosomeLength; ++i)
        {
            if(shiftChrBedRegionIdx)
            {
                final int currentRepeatStart = i - currentRepeatLength;
                // Note that chrBedRegions are indexed from one, not zero.
                while(chrBedRegionIdx < mChrBedRegions.size() && mChrBedRegions.get(chrBedRegionIdx).end() - 1 < currentRepeatStart)
                {
                    ++chrBedRegionIdx;
                }
                shiftChrBedRegionIdx = false;
            }

            final char base = Character.toUpperCase(refBases.charAt(i));
            if(history.size() <= 1 && base == 'N')
            {
                history.clear();
                currentRepeatLength = 0;
                shiftChrBedRegionIdx = true;
                continue;
            }

            if(history.isEmpty())
            {
                history.add(base);
                currentRepeatLength = 1;
                continue;
            }

            if(history.size() == 1)
            {
                history.add(base);
                currentRepeatLength = 2;
                continue;
            }

            if(history.get(0).equals(base))
            {
                history.set(0, history.get(1));
                history.set(1, base);
                ++currentRepeatLength;
                continue;
            }

            shiftChrBedRegionIdx = true;
            if(currentRepeatLength >= mConfig.MinRepeatBases)
            {
                final int repeatEnd = i - 1;
                final int repeatStart = repeatEnd - currentRepeatLength + 1;
                foundRepeat(repeatStart, repeatEnd, history, chrBedRegionIdx);
            }

            if(base == 'N')
            {
                history.clear();
                currentRepeatLength = 0;
                continue;
            }

            history.set(0, history.get(1));
            history.set(1, base);
            currentRepeatLength = 2;
        }

        if(currentRepeatLength >= mConfig.MinRepeatBases)
        {
            final int repeatEnd = chromosomeLength - 1;
            final int repeatStart = repeatEnd - currentRepeatLength + 1;
            foundRepeat(repeatStart, repeatEnd, history, chrBedRegionIdx);
        }
    }

    private void foundRepeat(final int repeatStart, final int repeatEnd, final List<Character> repeatBases, final int chrBedRegionIdx)
    {
        final ChrBaseRegion repeatRegion = new ChrBaseRegion(mChromosome, repeatStart + 1, repeatEnd + 1);
        RepeatInfo repeatInfo;
        if(mBedAnnotations)
        {
            final List<BaseRegion> containedWithinBedRegions = Lists.newArrayList();
            final List<BaseRegion> partialOverlapsWithBedRegions = Lists.newArrayList();
            for(int i = chrBedRegionIdx; i < mChrBedRegions.size() && mChrBedRegions.get(i).start() <= repeatRegion.end(); ++i)
            {
                final BaseRegion chrBedRegion = mChrBedRegions.get(i);
                if(chrBedRegion.containsRegion(repeatRegion))
                {
                    containedWithinBedRegions.add(chrBedRegion);
                }
                else if(chrBedRegion.overlaps(repeatRegion))
                {
                    partialOverlapsWithBedRegions.add(chrBedRegion);
                }
            }

            repeatInfo = new RepeatInfo(repeatBases, repeatRegion, containedWithinBedRegions, partialOverlapsWithBedRegions);
        }
        else
        {
            repeatInfo = new RepeatInfo(repeatBases, repeatRegion);
        }

        final boolean singleNuclRepeat = repeatBases.get(0).equals(repeatBases.get(1));
        if(singleNuclRepeat)
        {
            ++mSingleNuclRepeatCount;
            mTotalSingleNuclRepeatBases += repeatRegion.baseLength();
        }
        else
        {
            ++mDoubleNuclRepeatCount;
            mTotalDoubleNuclRepeatBases += repeatRegion.baseLength();
        }

        RepeatFinder.writeRepeat(mWriter, repeatInfo);
    }
}
