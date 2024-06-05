package com.hartwig.hmftools.esvee.depth;

import static java.lang.Math.abs;
import static java.lang.Math.max;
import static java.lang.Math.min;
import static java.lang.String.format;

import static com.hartwig.hmftools.common.genome.region.Orientation.FORWARD;
import static com.hartwig.hmftools.common.genome.region.Orientation.REVERSE;
import static com.hartwig.hmftools.common.region.BaseRegion.positionWithin;
import static com.hartwig.hmftools.common.region.BaseRegion.positionsOverlap;
import static com.hartwig.hmftools.common.bam.CigarUtils.leftSoftClipped;
import static com.hartwig.hmftools.common.bam.CigarUtils.rightSoftClipped;
import static com.hartwig.hmftools.common.sv.SvVcfTags.ALLELE_FRACTION;
import static com.hartwig.hmftools.common.sv.SvVcfTags.REF_DEPTH;
import static com.hartwig.hmftools.common.sv.SvVcfTags.REF_DEPTH_PAIR;
import static com.hartwig.hmftools.common.sv.SvVcfTags.TOTAL_FRAGS;
import static com.hartwig.hmftools.common.utils.PerformanceCounter.NANOS_IN_SECOND;
import static com.hartwig.hmftools.common.variant.CommonVcfTags.getGenotypeAttributeAsInt;
import static com.hartwig.hmftools.esvee.AssemblyConfig.SV_LOGGER;
import static com.hartwig.hmftools.esvee.prep.PrepConstants.DEFAULT_MAX_FRAGMENT_LENGTH;

import java.io.File;
import java.util.List;
import java.util.Map;
import java.util.StringJoiner;
import java.util.concurrent.Callable;

import com.google.common.annotations.VisibleForTesting;
import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import com.hartwig.hmftools.common.genome.region.Orientation;
import com.hartwig.hmftools.common.region.ChrBaseRegion;
import com.hartwig.hmftools.common.bam.BamSlicer;
import com.hartwig.hmftools.common.bam.SupplementaryReadData;
import com.hartwig.hmftools.common.utils.PerformanceCounter;

import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;
import htsjdk.variant.variantcontext.Genotype;
import htsjdk.variant.variantcontext.VariantContext;

public class DepthTask implements Callable
{
    private final DepthConfig mConfig;
    private final Map<String,Integer> mSampleVcfGenotypeIds;
    private final List<VariantContext> mVariantsList;
    private final List<VariantInfo> mVariantInfoList;
    private final String mChromosome;

    private final List<SamReader> mSamReaders;
    private final BamSlicer mBamSlicer;

    private final Map<String, ReadGroup> mReadGroups;

    private final SliceRegionState mSliceRegionState;
    private int mCurrentSampleIndex; // the sample / BAM being spliced

    private int mTotalReadCount;
    private int mCacheRecordCounter;
    private final PerformanceCounter mPerfCounter;

    public DepthTask(final String chromosome, final DepthConfig config, final Map<String,Integer> sampleVcfGenotypeIds)
    {
        mConfig = config;
        mChromosome = chromosome;
        mSampleVcfGenotypeIds = sampleVcfGenotypeIds;

        mVariantsList = Lists.newArrayList();
        mVariantInfoList = Lists.newArrayList();
        mTotalReadCount = 0;
        mCacheRecordCounter = 0;
        mSliceRegionState = new SliceRegionState();

        mReadGroups = Maps.newHashMap();
        mSamReaders = Lists.newArrayList();
        mBamSlicer = new BamSlicer(0, false, true, false);

        if(!mConfig.BamFiles.isEmpty())
        {
            for(String bamFile : mConfig.BamFiles)
            {
                mSamReaders.add(SamReaderFactory.makeDefault()
                        .validationStringency(mConfig.BamStringency)
                        .referenceSequence(new File(mConfig.RefGenome)).open(new File(bamFile)));
            }
        }

        mCurrentSampleIndex = 0;

        mPerfCounter = new PerformanceCounter("Slice");
    }

    public String chromosome() { return mChromosome; }

    public void addVariants(final List<VariantContext> variants)
    {
        List<Integer> genotypeIds = Lists.newArrayList();
        mConfig.Samples.forEach(x -> genotypeIds.add(mSampleVcfGenotypeIds.get(x)));

        for(VariantContext variant : variants)
        {
            mVariantsList.add(variant);
            mVariantInfoList.add(new VariantInfo(variant, genotypeIds, mConfig.VafCap));
        }
    }

    public List<VariantContext> variants() { return mVariantsList; }
    public PerformanceCounter getPerfCounter() { return mPerfCounter; }

    @Override
    public Long call()
    {
        SV_LOGGER.info("chr({}) processing {} variants", mChromosome, mVariantsList.size());

        // process the set of variants by grouping them into those with close positions where they may be able to share
        // the same reads from a wider slice
        int processed = 0;
        int index = 0;
        while(index < mVariantInfoList.size())
        {
            VariantInfo variant = mVariantInfoList.get(index);
            int posStart = variant.Position;

            mSliceRegionState.reset();
            mSliceRegionState.addVariant(variant);

            int posEnd = posStart;
            int nextIndex = index + 1;
            while(nextIndex < mVariantsList.size())
            {
                VariantInfo nextVariant = mVariantInfoList.get(nextIndex);
                if(nextVariant.Position - posEnd > mConfig.ProximityDistance)
                    break;

                posEnd = nextVariant.Position;
                mSliceRegionState.addVariant(nextVariant);
                ++nextIndex;
            }

            sliceSampleBams();

            index += mSliceRegionState.variantCount();

            processed += mSliceRegionState.variantCount();

            if((processed % 1000) == 0)
            {
                SV_LOGGER.debug("chr({}) processed {} variants", mChromosome, processed);
            }

            mCacheRecordCounter += mReadGroups.size();
            mReadGroups.clear();

            if(mCacheRecordCounter > READ_CACHE_CLEAR_COUNT)
            {
                // SV_LOGGER.debug("chr({}) read-group cache count({}) exceeds threshold", mChromosome, mCacheRecordCounter);
                mCacheRecordCounter = 0;
            }
        }

        // all current variants have had reads assigned to each sample, so now tally up their counts
        String refVcfTag = mConfig.getVcfTag(REF_DEPTH);
        String refPairVcfTag = mConfig.getVcfTag(REF_DEPTH_PAIR);

        for(int i = 0; i < mVariantsList.size(); ++i)
        {
            VariantContext variant = mVariantsList.get(i);
            VariantInfo variantInfo = mVariantInfoList.get(i);

            for(int s = 0; s < mConfig.Samples.size(); ++s)
            {
                String sampleId = mConfig.Samples.get(s);
                RefSupportCounts sampleCounts = variantInfo.SampleSupportCounts[s];
                int genotypeIndex = mSampleVcfGenotypeIds.get(sampleId);

                Genotype genotype = variant.getGenotype(genotypeIndex);

                if(genotype.getExtendedAttributes() == null || genotype.getExtendedAttributes().isEmpty())
                    continue;

                genotype.getExtendedAttributes().put(refVcfTag, sampleCounts.RefSupport);
                genotype.getExtendedAttributes().put(refPairVcfTag, sampleCounts.RefPairSupport);

                int variantFrags = getGenotypeAttributeAsInt(genotype, TOTAL_FRAGS, 0);

                double total = variantFrags + sampleCounts.total();
                double af = variantFrags / total;

                genotype.getExtendedAttributes().put(ALLELE_FRACTION, af);
            }

            RefSupportCounts totalCounts = variantInfo.totalSupport();
            setRefDepthValue(variant, totalCounts.RefSupport, refVcfTag);
            setRefDepthValue(variant, totalCounts.RefPairSupport, refPairVcfTag);
        }

        SV_LOGGER.info("chr({}) complete for {} variants, total reads({})", mChromosome, processed, mTotalReadCount);
        mReadGroups.clear();

        return (long)0;
    }

    private static final int READ_CACHE_CLEAR_COUNT = 100000;

    private void sliceSampleBams()
    {
        mPerfCounter.start();

        // retrieve the depth to set these 2 values in the VCF:
        // REF = reads aligned directly over the breakend (excluding reads that would support the junction)
        // REFPAIR = fragments with 1 read aligned to the left and 1 read aligned to the right of the breakend with proper orientation and
        // FragmentSize < max size from distribution

        // retrieve reads which may start and end before the first variant to capture read pair counts ie fragments which span the variants
        ChrBaseRegion region = new ChrBaseRegion(
                mChromosome,
                mSliceRegionState.PositionMin - DEFAULT_MAX_FRAGMENT_LENGTH,
                mSliceRegionState.PositionMax + DEFAULT_MAX_FRAGMENT_LENGTH);

        int readGroupTotal = 0;

        List<Double> times = Lists.newArrayList();
        List<Integer> readCounts = Lists.newArrayList();
        long startTime = 0;
        int readCount = 0;

        for(int i = 0; i < mSamReaders.size(); ++i)
        {
            mCurrentSampleIndex = i;

            SamReader samReader = mSamReaders.get(i);

            mReadGroups.clear();
            mSliceRegionState.resetUncappedVariants();

            startTime = System.nanoTime();
            readCount = mTotalReadCount;

            SV_LOGGER.trace("sample({}) slice for {} variants", mConfig.Samples.get(i), mSliceRegionState.variantCount());
            mBamSlicer.slice(samReader, region, this::processRead);

            times.add((System.nanoTime() - startTime)/NANOS_IN_SECOND);
            readCounts.add(mTotalReadCount - readCount);

            mReadGroups.values().forEach(x -> processReadGroup(x));
        }

        mPerfCounter.stop();

        if(mConfig.PerfLogTime > 0 &&  mPerfCounter.getLastTime() > mConfig.PerfLogTime)
        {
            StringJoiner sjTimes = new StringJoiner(",");
            StringJoiner sjCounts = new StringJoiner(",");
            times.forEach(x -> sjTimes.add(format("%.3f", x)));
            readCounts.forEach(x -> sjCounts.add(format("%d", x)));
            SV_LOGGER.debug("chr({}) slice({}) high depth retrieval time({}) totalFrags({}) times({}) readCounts({})",
                    mChromosome, mSliceRegionState, format("%.3f", mPerfCounter.getLastTime()), readGroupTotal,
                    sjTimes.toString(), sjCounts.toString());
        }
    }

    private void processRead(final SAMRecord read)
    {
        ++mTotalReadCount;

        int maxSlicePosition = mSliceRegionState.PositionMax + DEFAULT_MAX_FRAGMENT_LENGTH;

        ReadGroup readGroup = mReadGroups.get(read.getReadName());

        if(readGroup != null)
        {
            readGroup.Reads.add(read);

            if(!readGroup.WaitForAll)
            {
                if(read.getSupplementaryAlignmentFlag())
                {
                    readGroup.WaitForAll = true;
                }
                else
                {
                    // check again for a supplementary
                    SupplementaryReadData suppReadData = SupplementaryReadData.extractAlignment(read);

                    if(suppReadData != null && suppReadData.Chromosome.equals(mChromosome)
                            && positionWithin(suppReadData.Position, mSliceRegionState.PositionMin, maxSlicePosition))
                    {
                        readGroup.WaitForAll = true;
                    }
                }
            }

            if(!readGroup.WaitForAll)
            {
                processReadGroup(readGroup);
                mReadGroups.remove(read.getReadName());
            }

            return;
        }

        // check for any support for this read against the current variants
        boolean isRelevant = false;
        boolean anyUncapped = false;

        int index = 0;
        while(index < mSliceRegionState.UncappedVariants.size())
        {
            VariantInfo variantInfo = mSliceRegionState.UncappedVariants.get(index);

            RefSupportCounts variantSampleCounts = variantInfo.SampleSupportCounts[mCurrentSampleIndex];

            if(variantSampleCounts.exceedsMaxDepth())
            {
                mSliceRegionState.UncappedVariants.remove(index);
                mSliceRegionState.MinPositionIndex = max(mSliceRegionState.MinPositionIndex - 1, 0);
                continue;
            }

            anyUncapped = true;

            if(isRelevantRead(read, variantInfo))
            {
                isRelevant = true;
                break;
            }

            ++index;
        }

        // stop slicing this region for this sample if all variants have reached their VAF cap
        if(!anyUncapped)
        {
            SV_LOGGER.debug("chr({}) variants({}) range({} - {}) sampleIndex({}) VAF cap exceeded for all variants",
                    mChromosome, mSliceRegionState.variantCount(), mSliceRegionState.PositionMin, mSliceRegionState.PositionMax,
                    mCurrentSampleIndex);

            mBamSlicer.haltProcessing();
            return;
        }

        if(!isRelevant)
            return;

        // determine if mate or supp reads are expected within this current slice region
        boolean expectSupplementaries = read.getSupplementaryAlignmentFlag();
        boolean expectMate = false;

        if(read.getReadPairedFlag() && !read.getMateUnmappedFlag()
        && read.getMateReferenceIndex() == read.getReferenceIndex() && read.getMateAlignmentStart() <= maxSlicePosition)
        {
            expectMate = true;
        }

        SupplementaryReadData suppReadData = SupplementaryReadData.extractAlignment(read);

        if(suppReadData != null && suppReadData.Chromosome.equals(mChromosome)
        && positionWithin(suppReadData.Position, mSliceRegionState.PositionMin, maxSlicePosition))
        {
            expectSupplementaries = true;
        }

        if(!expectMate && !expectSupplementaries)
        {
            processSingleRead(read);
        }
        else
        {
            mReadGroups.put(read.getReadName(), new ReadGroup(read, expectSupplementaries));
        }
    }

    private boolean isRelevantRead(final SAMRecord read, VariantInfo variantInfo)
    {
        // ignore reads which cannot span the current variant(s)
        if(read.getAlignmentEnd() < variantInfo.PositionMin)
        {
            if(!read.getReadPairedFlag() || read.getMateUnmappedFlag() || read.getMateReferenceIndex() != read.getReferenceIndex())
                return false;

            // both reads before the variants
            if(read.getMateAlignmentStart() + read.getReadBases().length < variantInfo.PositionMin)
                return false;
        }
        else if(read.getAlignmentStart() > variantInfo.PositionMax)
        {
            if(!read.getReadPairedFlag() || read.getMateUnmappedFlag() || read.getMateReferenceIndex() != read.getReferenceIndex())
                return false;

            return false;
        }

        return true;
    }

    private static final int READ_POSITION_MARGIN = 100;

    private void processSingleRead(final SAMRecord read)
    {
        if(mSliceRegionState.UncappedVariants.size() < 10)
        {
            processReadGroup(new ReadGroup(read, false));
            return;
        }

        int newMinPositionIndex = mSliceRegionState.MinPositionIndex;

        for(int i = mSliceRegionState.MinPositionIndex; i < mSliceRegionState.UncappedVariants.size(); ++i)
        {
            VariantInfo variantInfo = mSliceRegionState.UncappedVariants.get(i);
            if(read.getAlignmentStart() > variantInfo.PositionMax + READ_POSITION_MARGIN)
            {
                newMinPositionIndex = i + 1;
                continue;
            }
            else if(read.getAlignmentEnd() < variantInfo.PositionMin - READ_POSITION_MARGIN)
            {
                break;
            }

            RefSupportCounts variantSampleCounts = variantInfo.SampleSupportCounts[mCurrentSampleIndex];
            checkReadGroupSupport(variantInfo, variantSampleCounts, new ReadGroup(read, false));
        }

        // record the starting index for the set of current variants
        mSliceRegionState.MinPositionIndex = newMinPositionIndex;
    }

    private void processReadGroup(final ReadGroup readGroup)
    {
        // find the variants that this group overlaps
        int groupMinPosition = readGroup.Reads.stream().mapToInt(x -> x.getAlignmentStart()).min().orElse(0);
        int groupMaxPosition = readGroup.Reads.stream().mapToInt(x -> x.getAlignmentEnd()).max().orElse(0);

        for(VariantInfo variantInfo : mSliceRegionState.UncappedVariants)
        {
            if(groupMinPosition > variantInfo.PositionMax + READ_POSITION_MARGIN)
                continue;
            else if(groupMaxPosition < variantInfo.PositionMin - READ_POSITION_MARGIN)
                break;

            RefSupportCounts variantSampleCounts = variantInfo.SampleSupportCounts[mCurrentSampleIndex];
            checkReadGroupSupport(variantInfo, variantSampleCounts, readGroup);
        }
    }

    private void checkReadGroupSupport(final VariantInfo variant, RefSupportCounts supportCounts, final ReadGroup readGroup)
    {
        boolean readSupportsRef = false;
        boolean hasLowerPosRead = false;
        boolean hasUpperPosRead = false;
        int strandCount = 0;
        boolean matchesJunction = false;
        int readGroupPosMin = 0;
        int readGroupPosMax = 0;

        for(SAMRecord read : readGroup.Reads)
        {
            int readStart = read.getAlignmentStart();
            int readEnd = read.getAlignmentEnd();
            boolean isSupplementary = read.getSupplementaryAlignmentFlag();

            if(!isSupplementary)
            {
                readGroupPosMin = readGroupPosMin == 0 ? readStart : min(readStart, readGroupPosMin);
                readGroupPosMax = max(readEnd, readGroupPosMax);
            }

            // check for an exact SC match
            if((variant.Orient.isReverse() && positionWithin(readStart, variant.PositionMin, variant.PositionMax) && leftSoftClipped(read))
            || (variant.Orient.isForward() && positionWithin(readEnd, variant.PositionMin, variant.PositionMax)) && rightSoftClipped(read))
            {
                SV_LOGGER.trace("var({}) pos({}-{}) read({}-{}) id({}) at junction",
                        variant.Position, variant.PositionMin, variant.PositionMax, readStart, readEnd, read.getReadName());
                matchesJunction = true;
                break;
            }

            if(!isSupplementary && readGroup.Reads.size() > 1)
            {
                Orientation orientation = !read.getReadNegativeStrandFlag() ? FORWARD : REVERSE;

                if(orientation.isForward() && readEnd <= max(variant.Position, variant.PositionMax) && !hasLowerPosRead
                && abs(read.getInferredInsertSize()) < DEFAULT_MAX_FRAGMENT_LENGTH)
                {
                    hasLowerPosRead = true;
                    strandCount += read.getReadNegativeStrandFlag() ? -1 : 1;
                }
                else if(orientation.isReverse() && readStart >= min(variant.Position, variant.PositionMin) && !hasUpperPosRead
                && abs(read.getInferredInsertSize()) < DEFAULT_MAX_FRAGMENT_LENGTH)
                {
                    hasUpperPosRead = true;
                    strandCount += read.getReadNegativeStrandFlag() ? -1 : 1;
                }
            }

            if(positionsOverlap(variant.PositionMin, variant.PositionMax, readStart, readEnd))
            {
                //SV_LOGGER.trace("var({}) pos({}-{}) read({}-{}) id({}) has ref support",
                //        variant.Position, variant.PositionMin, variant.PositionMax, readStart, readEnd, read.getReadName());
                readSupportsRef = true;
            }
        }

        if(matchesJunction)
            return;

        if(readSupportsRef)
        {
            ++supportCounts.RefSupport;
        }
        else if(hasLowerPosRead && hasUpperPosRead && strandCount == 0)
        {
            ++supportCounts.RefPairSupport;

            SV_LOGGER.trace("var({}) pos({}-{}) fragment(id={} {}-{}) has ref-pair support",
                    variant.Position, variant.PositionMin, variant.PositionMax, readGroup.id(), readGroupPosMin, readGroupPosMax);
        }

        if(supportCounts.exceedsMaxDepth())
        {
            SV_LOGGER.trace("var({}:{}) sampleIndex({}) ref limit({}) reached with support(ref={} pair={})",
                    mChromosome, variant.Position, mCurrentSampleIndex,
                    supportCounts.VafCap, supportCounts.RefSupport, supportCounts.RefPairSupport);
        }
    }

    private void setRefDepthValue(final VariantContext variant, int refCount, final String vcfTag)
    {
        if(variant.hasAttribute(vcfTag))
            variant.getCommonInfo().removeAttribute(vcfTag);

        variant.getCommonInfo().putAttribute(vcfTag, refCount);
    }

    @VisibleForTesting
    public void reset()
    {
        mVariantInfoList.clear();
        mVariantsList.clear();
        mReadGroups.clear();
        mCurrentSampleIndex = 0;
        mCacheRecordCounter = 0;
        mTotalReadCount = 0;
        mSliceRegionState.reset();
        mSliceRegionState.reset();
    }

    @VisibleForTesting
    public void processSamRecord(final SAMRecord read)
    {
        processRead(read);
    }

    @VisibleForTesting
    public Map<String,ReadGroup> readGroups() { return mReadGroups; }

    @VisibleForTesting
    public List<VariantInfo> variantInfos() { return mVariantInfoList; }

    @VisibleForTesting
    public void addSliceVariants(final List<VariantInfo> variants)
    {
        mSliceRegionState.reset();
        variants.forEach(x -> mSliceRegionState.addVariant(x));
        mSliceRegionState.resetUncappedVariants();
    }

    @VisibleForTesting
    public SliceRegionState sliceRegionState() { return mSliceRegionState; }
    }
