package com.hartwig.hmftools.svprep.depth;

import static java.lang.Math.abs;
import static java.lang.Math.max;
import static java.lang.Math.min;
import static java.lang.String.format;

import static com.hartwig.hmftools.common.samtools.CigarUtils.leftSoftClipped;
import static com.hartwig.hmftools.common.samtools.CigarUtils.rightSoftClipped;
import static com.hartwig.hmftools.common.samtools.SamRecordUtils.SUPPLEMENTARY_ATTRIBUTE;
import static com.hartwig.hmftools.common.sv.StructuralVariantFactory.ALLELE_FRACTION;
import static com.hartwig.hmftools.common.sv.StructuralVariantFactory.REFERENCE_BREAKEND_READPAIR_COVERAGE;
import static com.hartwig.hmftools.common.sv.StructuralVariantFactory.REFERENCE_BREAKEND_READ_COVERAGE;
import static com.hartwig.hmftools.common.sv.StructuralVariantFactory.VARIANT_FRAGMENT_BREAKEND_COVERAGE;
import static com.hartwig.hmftools.common.sv.StructuralVariantFactory.VARIANT_FRAGMENT_BREAKPOINT_COVERAGE;
import static com.hartwig.hmftools.common.utils.PerformanceCounter.NANOS_IN_SECOND;
import static com.hartwig.hmftools.common.utils.sv.BaseRegion.positionWithin;
import static com.hartwig.hmftools.common.utils.sv.BaseRegion.positionsOverlap;
import static com.hartwig.hmftools.common.utils.sv.SvCommonUtils.NEG_ORIENT;
import static com.hartwig.hmftools.common.utils.sv.SvCommonUtils.POS_ORIENT;
import static com.hartwig.hmftools.common.variant.CommonVcfTags.getGenotypeAttributeAsInt;
import static com.hartwig.hmftools.svprep.SvCommon.SV_LOGGER;
import static com.hartwig.hmftools.svprep.SvConstants.DEFAULT_MAX_FRAGMENT_LENGTH;
import static com.hartwig.hmftools.svprep.depth.DepthConfig.MAX_GAP;

import java.io.File;
import java.util.List;
import java.util.Map;
import java.util.StringJoiner;
import java.util.concurrent.Callable;

import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import com.hartwig.hmftools.common.samtools.BamSlicer;
import com.hartwig.hmftools.common.samtools.SupplementaryReadData;
import com.hartwig.hmftools.common.utils.PerformanceCounter;
import com.hartwig.hmftools.common.utils.sv.ChrBaseRegion;

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

        mBamSlicer = new BamSlicer(0, false, true, false);
        mReadGroups = Maps.newHashMap();

        mSamReaders = Lists.newArrayList();

        for(String bamFile : mConfig.BamFiles)
        {
            mSamReaders.add(SamReaderFactory.makeDefault().referenceSequence(new File(mConfig.RefGenome)).open(new File(bamFile)));
        }

        mCurrentSampleIndex = 0;

        mPerfCounter = new PerformanceCounter("Slice");
    }

    public String chromosome() { return mChromosome; }

    public void addVariants(final List<VariantContext> variants)
    {
        for(VariantContext variant : variants)
        {
            mVariantsList.add(variant);
            mVariantInfoList.add(new VariantInfo(variant, mConfig.Samples.size(), mConfig.VafCap));
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
            mSliceRegionState.addVariant(index, variant);

            int posEnd = posStart;
            int nextIndex = index + 1;
            while(nextIndex < mVariantsList.size())
            {
                VariantInfo nextVariant = mVariantInfoList.get(nextIndex);
                if(nextVariant.Position - posEnd > MAX_GAP)
                    break;

                posEnd = nextVariant.Position;
                mSliceRegionState.addVariant(nextIndex, nextVariant);
                ++nextIndex;
            }

            sliceSampleBams();

            index += mSliceRegionState.VariantCount;

            processed += mSliceRegionState.VariantCount;

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
                System.gc();
            }
        }

        // all current variants have had reads assigned to each sample, so now tally up their counts
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

                genotype.getExtendedAttributes().put(REFERENCE_BREAKEND_READ_COVERAGE, sampleCounts.RefSupport);
                genotype.getExtendedAttributes().put(REFERENCE_BREAKEND_READPAIR_COVERAGE, sampleCounts.RefPairSupport);

                int variantFrags = getGenotypeAttributeAsInt(genotype, VARIANT_FRAGMENT_BREAKPOINT_COVERAGE, 0) +
                        getGenotypeAttributeAsInt(genotype, VARIANT_FRAGMENT_BREAKEND_COVERAGE, 0);

                double total = variantFrags + sampleCounts.total();
                double af = variantFrags / total;

                genotype.getExtendedAttributes().put(ALLELE_FRACTION, af);
            }

            RefSupportCounts totalCounts = variantInfo.totalSupport();
            setRefDepthValue(variant, totalCounts.RefSupport, REFERENCE_BREAKEND_READ_COVERAGE);
            setRefDepthValue(variant, totalCounts.RefPairSupport, REFERENCE_BREAKEND_READPAIR_COVERAGE);
        }

        SV_LOGGER.info("chr({}) complete for {} variants, total reads({})", mChromosome, processed, mTotalReadCount);
        mReadGroups.clear();
        System.gc();

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

            startTime = System.nanoTime();
            readCount = mTotalReadCount;

            mBamSlicer.slice(samReader, Lists.newArrayList(region), this::processSamRecord);

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

    private void processSamRecord(final SAMRecord read)
    {
        ++mTotalReadCount;

        ReadGroup readGroup = mReadGroups.get(read.getReadName());

        if(readGroup != null)
        {
            readGroup.Reads.add(read);

            if(!readGroup.WaitForAll)
            {
                processReadGroup(readGroup);
                mReadGroups.remove(read.getReadName());
            }

            return;
        }

        // check for any support for this read against the current variants
        boolean isRelevant = false;

        for(int i = mSliceRegionState.VariantIndexStart; i <= mSliceRegionState.VariantIndexEnd; ++i)
        {
            if(isRelevantRead(read, mVariantInfoList.get(i)))
            {
                isRelevant = true;
                break;
            }
        }

        if(!isRelevant)
            return;

        // determine if mate or supp reads are expected within this current slice region
        boolean expectSupplementaries = false;
        boolean expectMate = false;

        if(read.getReadPairedFlag() && !read.getMateUnmappedFlag()
        && read.getMateReferenceIndex() == read.getReferenceIndex() && read.getMateAlignmentStart() <= mSliceRegionState.PositionMax)
        {
            expectMate = true;
        }

        SupplementaryReadData suppReadData = SupplementaryReadData.from(read);

        if(suppReadData != null && suppReadData.Chromosome.equals(mChromosome)
        && positionWithin(suppReadData.Position, mSliceRegionState.PositionMin, mSliceRegionState.PositionMax))
        {
            expectSupplementaries = true;
        }

        if(!expectMate && !expectSupplementaries)
        {
            processReadGroup(new ReadGroup(read, false));
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
            if(read.getMateUnmappedFlag() || read.getMateReferenceIndex() != read.getReferenceIndex())
                return false;

            // both reads before the variants
            if(read.getMateAlignmentStart() + read.getReadBases().length < variantInfo.PositionMin)
                return false;
        }
        else if(read.getAlignmentStart() > variantInfo.PositionMax)
        {
            if(read.getMateUnmappedFlag() || read.getMateReferenceIndex() != read.getReferenceIndex())
                return false;

            // the mate's position isn't checked if it's an earlier read
            //if(read.getMateAlignmentStart() > variantInfo.PositionMax)

            return false;
        }

        return true;
    }

    private void processReadGroup(final ReadGroup readGroup)
    {
        // find the variants that this group overlaps
        for(int i = mSliceRegionState.VariantIndexStart; i <= mSliceRegionState.VariantIndexEnd; ++i)
        {
            VariantInfo variantInfo = mVariantInfoList.get(i);
            RefSupportCounts variantSampleCounts = variantInfo.SampleSupportCounts[mCurrentSampleIndex];

            if(variantInfo.RefFragsCap > 0 && variantSampleCounts.total() > variantInfo.RefFragsCap)
                continue;

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
            if((variant.Orientation == NEG_ORIENT && positionWithin(readStart, variant.PositionMin, variant.PositionMax) && leftSoftClipped(read))
            || (variant.Orientation == POS_ORIENT && positionWithin(readEnd, variant.PositionMin, variant.PositionMax)) && rightSoftClipped(read))
            {
                SV_LOGGER.trace("var({}) pos({}-{}) read({}-{}) id({}) at junction",
                        variant.Position, variant.PositionMin, variant.PositionMax, readStart, readEnd, read.getReadName());
                matchesJunction = true;
                break;
            }

            if(!isSupplementary)
            {
                byte orientation = !read.getReadNegativeStrandFlag() ? POS_ORIENT : NEG_ORIENT;

                if(orientation == POS_ORIENT && readEnd <= max(variant.Position, variant.PositionMax) && !hasLowerPosRead
                        && abs(read.getInferredInsertSize()) < DEFAULT_MAX_FRAGMENT_LENGTH)
                {
                    hasLowerPosRead = true;
                    strandCount += read.getReadNegativeStrandFlag() ? -1 : 1;
                }
                else if(orientation == NEG_ORIENT && readStart >= min(variant.Position, variant.PositionMin) && !hasUpperPosRead
                        && abs(read.getInferredInsertSize()) < DEFAULT_MAX_FRAGMENT_LENGTH)
                {
                    hasUpperPosRead = true;
                    strandCount += read.getReadNegativeStrandFlag() ? -1 : 1;
                }
            }

            if(positionsOverlap(variant.PositionMin, variant.PositionMax, readStart, readEnd))
            {
                SV_LOGGER.trace("var({}) pos({}-{}) read({}-{}) id({}) has ref support",
                        variant.Position, variant.PositionMin, variant.PositionMax, readStart, readEnd, read.getReadName());
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

        if(variant.RefFragsCap > 0 && supportCounts.total() > variant.RefFragsCap)
        {
            SV_LOGGER.debug("var({}:{}) varFrags({}) ref limit({}) reached with support(ref={} pair={})",
                    mChromosome, variant.Position, variant.FragmentCount, variant.RefFragsCap,
                    supportCounts.RefSupport, supportCounts.RefPairSupport);
        }
    }

    private void setRefDepthValue(final VariantContext variant, int refCount, final String vcfTag)
    {
        if(variant.hasAttribute(vcfTag))
            variant.getCommonInfo().removeAttribute(vcfTag);

        variant.getCommonInfo().putAttribute(vcfTag, refCount);
    }
}
