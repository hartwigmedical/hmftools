package com.hartwig.hmftools.svprep.depth;

import static java.lang.Math.abs;
import static java.lang.Math.max;
import static java.lang.Math.min;
import static java.lang.String.format;

import static com.hartwig.hmftools.common.sv.StructuralVariantFactory.CIPOS;
import static com.hartwig.hmftools.common.sv.StructuralVariantFactory.REFERENCE_BREAKEND_READPAIR_COVERAGE;
import static com.hartwig.hmftools.common.sv.StructuralVariantFactory.REFERENCE_BREAKEND_READ_COVERAGE;
import static com.hartwig.hmftools.common.sv.StructuralVariantFactory.VARIANT_FRAGMENT_BREAKEND_COVERAGE;
import static com.hartwig.hmftools.common.sv.StructuralVariantFactory.VARIANT_FRAGMENT_BREAKPOINT_COVERAGE;
import static com.hartwig.hmftools.common.sv.StructuralVariantFactory.isSingleBreakend;
import static com.hartwig.hmftools.common.utils.sv.BaseRegion.positionWithin;
import static com.hartwig.hmftools.common.utils.sv.BaseRegion.positionsOverlap;
import static com.hartwig.hmftools.common.utils.sv.SvCommonUtils.NEG_ORIENT;
import static com.hartwig.hmftools.common.utils.sv.SvCommonUtils.POS_ORIENT;
import static com.hartwig.hmftools.svprep.SvCommon.SV_LOGGER;
import static com.hartwig.hmftools.svprep.SvConstants.DEFAULT_MAX_FRAGMENT_LENGTH;
import static com.hartwig.hmftools.svprep.depth.DepthAnnotator.VCF_TAG_REFPAIR_GRIDSS;
import static com.hartwig.hmftools.svprep.depth.DepthAnnotator.VCF_TAG_REF_GRIDSS;
import static com.hartwig.hmftools.svprep.depth.DepthConfig.MAX_GAP;

import java.io.File;
import java.util.List;
import java.util.Map;
import java.util.concurrent.Callable;

import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import com.hartwig.hmftools.common.samtools.BamSlicer;
import com.hartwig.hmftools.common.utils.PerformanceCounter;
import com.hartwig.hmftools.common.utils.sv.ChrBaseRegion;

import org.apache.logging.log4j.Level;

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
    private final String mChromosome;
    private int mBamRecords;

    private final List<SamReader> mSamReaders;
    private final BamSlicer mBamSlicer;

    private final Map<String,SamReadGroup> mReadGroups;

    private final PerformanceCounter mPerfCounter;

    public DepthTask(final String chromosome, final DepthConfig config, final Map<String,Integer> sampleVcfGenotypeIds)
    {
        mConfig = config;
        mChromosome = chromosome;
        mSampleVcfGenotypeIds = sampleVcfGenotypeIds;

        mVariantsList = Lists.newArrayList();
        mBamRecords = 0;

        mBamSlicer = new BamSlicer(0, false, true, false);
        mReadGroups = Maps.newHashMap();

        mSamReaders = Lists.newArrayList();

        for(String bamFile : mConfig.BamFiles)
        {
            mSamReaders.add(SamReaderFactory.makeDefault().referenceSequence(new File(mConfig.RefGenome)).open(new File(bamFile)));
        }

        mPerfCounter = new PerformanceCounter("Slice");
    }

    public String chromosome() { return mChromosome; }
    public void addVariants(final List<VariantContext> variants)
    {
        mVariantsList.addAll(variants);
    }

    public List<VariantContext> variants() { return mVariantsList; }
    public PerformanceCounter getPerfCounter() { return mPerfCounter; }

    @Override
    public Long call()
    {
        SV_LOGGER.info("chr({}) processing {} variants", mChromosome, mVariantsList.size());

        int processed = 0;
        int index = 0;
        while(index < mVariantsList.size())
        {
            VariantContext variant = mVariantsList.get(index);
            int posStart = variant.getStart();

            List<VariantContext> variants = Lists.newArrayList(variant);

            int posEnd = posStart;
            int nextIndex = index + 1;
            while(nextIndex < mVariantsList.size())
            {
                VariantContext nextVariant = mVariantsList.get(nextIndex);
                if(nextVariant.getStart() - posEnd > MAX_GAP)
                    break;

                posEnd = nextVariant.getStart();
                variants.add(nextVariant);
                ++nextIndex;
            }

            retrieveDepth(variants, posStart, posEnd);
            index += variants.size();

            processed += variants.size();

            if((processed % 1000) == 0)
            {
                SV_LOGGER.debug("chr({}) processed {} variants", mChromosome, processed);
            }
        }

        SV_LOGGER.info("chr({}) complete for {} variants, total reads({})", mChromosome, processed, mBamRecords);
        mReadGroups.clear();
        System.gc();

        return (long)0;
    }

    private void retrieveDepth(final List<VariantContext> variants, int posStart, int posEnd)
    {
        mPerfCounter.start();

        // retrieve the depth to set these 2 values in the VCF:
        // REF = reads aligned directly over the breakend (excluding reads that would support the junction)
        // REFPAIR = fragments with 1 read aligned to the left and 1 read aligned to the right of the breakend with proper orientation and
        // FragmentSize < max size from distribution

        ChrBaseRegion region = new ChrBaseRegion(
                mChromosome, posStart - DEFAULT_MAX_FRAGMENT_LENGTH, posEnd + DEFAULT_MAX_FRAGMENT_LENGTH);

        List<RefSupportCounts> sampleTotalCounts = Lists.newArrayList();

        for(int i = 0; i < variants.size(); ++i)
        {
            sampleTotalCounts.add(new RefSupportCounts());
        }

        int readGroupTotal = 0;

        for(int i = 0; i < mConfig.Samples.size(); ++i)
        {
            String sampleId = mConfig.Samples.get(i);
            SamReader samReader = mSamReaders.get(i);

            mReadGroups.clear();

            mBamSlicer.slice(samReader, Lists.newArrayList(region), this::processSamRecord);

            readGroupTotal += mReadGroups.size();

            for(int j = 0; j < variants.size(); ++j)
            {
                VariantContext variant = variants.get(j);
                RefSupportCounts totalCounts = sampleTotalCounts.get(j);
                calculateVariantSupport(variant, sampleId, totalCounts);
            }
        }

        for(int j = 0; j < variants.size(); ++j)
        {
            RefSupportCounts totalCounts = sampleTotalCounts.get(j);
            VariantContext variant = variants.get(j);

            setRefDepthValue(variant, totalCounts.RefSupport, REFERENCE_BREAKEND_READ_COVERAGE, VCF_TAG_REF_GRIDSS);
            setRefDepthValue(variant, totalCounts.RefPairSupport, REFERENCE_BREAKEND_READPAIR_COVERAGE, VCF_TAG_REFPAIR_GRIDSS);
        }

        mPerfCounter.stop();

        if(mPerfCounter.getLastTime() > 2)
        {
            SV_LOGGER.warn("chr({}) span({}-{}) variants({}) high depth retrieval time({}) totalFrags({})",
                    mChromosome, posStart, posEnd, variants.size(), format("%.3f", mPerfCounter.getLastTime()), readGroupTotal);
        }
    }

    private void calculateVariantSupport(
            final VariantContext variant, final String sampleId, final RefSupportCounts totalCounts)
    {
        int variantPosition = variant.getStart();
        byte variantOrientation = getOrientation(variant);

        final int[] homology = {0, 0};

        if(variant.hasAttribute(CIPOS))
        {
            final List<Integer> ihompos = variant.getAttributeAsIntList(CIPOS, 0);
            homology[0] = ihompos.get(0);
            homology[1] = ihompos.get(1);
        }

        int varPosStart = variantPosition + homology[0];
        int varPosEnd = variantPosition + homology[1];

        int variantFragCount = max(
                variant.getAttributeAsInt(VARIANT_FRAGMENT_BREAKPOINT_COVERAGE, 0),
                variant.getAttributeAsInt(VARIANT_FRAGMENT_BREAKEND_COVERAGE, 0));

        RefSupportCounts sampleCounts = calculateSupport(variantPosition, varPosStart, varPosEnd, variantOrientation, variantFragCount);

        int vcfSampleIndex = mSampleVcfGenotypeIds.get(sampleId);
        Genotype genotype = variant.getGenotype(vcfSampleIndex);

        if(mConfig.LogDiffs || SV_LOGGER.isTraceEnabled())
        {
            int vcfRef = getGenotypeAttributeAsInt(genotype, REFERENCE_BREAKEND_READ_COVERAGE, 0);
            int vcfRefPair = getGenotypeAttributeAsInt(genotype, REFERENCE_BREAKEND_READPAIR_COVERAGE, 0);

            boolean hasDiffs = hasDiff(vcfRef, sampleCounts.RefSupport) || hasDiff(vcfRefPair, sampleCounts.RefPairSupport);

            Level level = hasDiffs ? Level.DEBUG : Level.TRACE;

            SV_LOGGER.log(level, "sample({}) var({}:{}) {} REF(vcf={} calc={}) REFPAIR(vcf={} calc={})",
                    sampleId, variant.getContig(), variantPosition, hasDiffs ? "has diffs" : "equal",
                    vcfRef, sampleCounts.RefSupport, vcfRefPair, sampleCounts.RefPairSupport);
        }

        totalCounts.RefSupport += sampleCounts.RefSupport;
        totalCounts.RefPairSupport += sampleCounts.RefPairSupport;

        setRefDepthValue(genotype, sampleCounts.RefSupport, REFERENCE_BREAKEND_READ_COVERAGE, VCF_TAG_REF_GRIDSS);
        setRefDepthValue(genotype, sampleCounts.RefPairSupport, REFERENCE_BREAKEND_READPAIR_COVERAGE, VCF_TAG_REFPAIR_GRIDSS);
    }

    private void setRefDepthValue(final Genotype genotype, int refCount, final String vcfTag, final String oldValueTag)
    {
        if(genotype.hasExtendedAttribute(vcfTag))
        {
            int oldValue = getGenotypeAttributeAsInt(genotype, vcfTag, 0);

            if(mConfig.WriteGridssRefValues)
                genotype.getExtendedAttributes().put(oldValueTag, oldValue);
        }

        genotype.getExtendedAttributes().put(vcfTag, refCount);
    }

    private void setRefDepthValue(final VariantContext variant, int refCount, final String vcfTag, final String oldValueTag)
    {
        if(variant.hasAttribute(vcfTag))
        {
            int oldValue = variant.getAttributeAsInt(vcfTag, 0);
            variant.getCommonInfo().removeAttribute(vcfTag);

            if(mConfig.WriteGridssRefValues)
                variant.getCommonInfo().putAttribute(oldValueTag, oldValue);
        }

        variant.getCommonInfo().putAttribute(vcfTag, refCount);
    }

    private static final int ABS_DIFF_MAX = 3;
    private static final double ABS_DIFF_PERC = 0.1;

    private boolean hasDiff(int value1, int value2)
    {
        int diff = abs(value1 - value2);
        double diffPerc = diff / (double)max(value1, value2);
        return diffPerc > ABS_DIFF_PERC && diff > ABS_DIFF_MAX;
    }

    private static int getGenotypeAttributeAsInt(final Genotype genotype, final String attribute, int defaultVaue)
    {
        Object value = genotype.getExtendedAttribute(attribute);
        return value == null ? defaultVaue : Integer.parseInt(value.toString());
    }

    private byte getOrientation(final VariantContext variant)
    {
        String alt = variant.getAlternateAllele(0).getDisplayString();

        if(isSingleBreakend(variant))
        {
            return alt.startsWith(".") ? NEG_ORIENT : POS_ORIENT;
        }
        else
        {
            return alt.startsWith("]") || alt.startsWith("[") ? NEG_ORIENT : POS_ORIENT;
        }
    }

    private void processSamRecord(final SAMRecord record)
    {
        ++mBamRecords;

        SamReadGroup readGroup = mReadGroups.get(record.getReadName());

        if(readGroup == null)
        {
            mReadGroups.put(record.getReadName(), new SamReadGroup(record));
            return;
        }

        readGroup.Reads.add(record);
    }

    private class RefSupportCounts
    {
        public int RefSupport = 0;
        public int RefPairSupport = 0;

        public RefSupportCounts() {}

        public int total() { return RefSupport + RefPairSupport; }
    }

    private RefSupportCounts calculateSupport(int variantPosition, int varPosMin, int varPosMax, byte variantOrientation, int varFrags)
    {
        RefSupportCounts counts = new RefSupportCounts();

        int refFragsCap = mConfig.VafCap > 0 ? (int)(varFrags / mConfig.VafCap) : 0;

        for(SamReadGroup readGroup : mReadGroups.values())
        {
            boolean readSupportsRef = false;
            boolean hasLowerPosRead = false;
            boolean hasUpperPosRead = false;
            int strandCount = 0;
            boolean matchesJunction = false;
            int readGroupPosMin = 0;
            int readGroupPosMax = 0;

            for(SAMRecord record : readGroup.Reads)
            {
                int readStart = record.getAlignmentStart();
                int readEnd = record.getAlignmentEnd();

                readGroupPosMin = readGroupPosMin == 0 ? readStart : min(readStart, readGroupPosMin);
                readGroupPosMax = max(readEnd, readGroupPosMax);

                // check for an exact SC match
                if((variantOrientation == NEG_ORIENT && positionWithin(readStart, varPosMin, varPosMax) && record.getCigar().isLeftClipped())
                || (variantOrientation == POS_ORIENT && positionWithin(readEnd, varPosMin, varPosMax)) && record.getCigar().isRightClipped())
                {
                    SV_LOGGER.trace("var({}) pos({}-{}) read({}-{}) id({}) at junction",
                            variantPosition, varPosMin, varPosMax, readStart, readEnd, record.getReadName());
                    matchesJunction = true;
                    break;
                }

                byte orientation = !record.getReadNegativeStrandFlag() ? POS_ORIENT : NEG_ORIENT;

                if(orientation == POS_ORIENT && readEnd <= max(variantPosition, varPosMax) && !hasLowerPosRead
                && abs(record.getInferredInsertSize()) < DEFAULT_MAX_FRAGMENT_LENGTH)
                {
                    hasLowerPosRead = true;
                    strandCount += record.getReadNegativeStrandFlag() ? -1 : 1;
                }
                else if(orientation == NEG_ORIENT && readStart >= min(variantPosition, varPosMin) && !hasUpperPosRead
                && abs(record.getInferredInsertSize()) < DEFAULT_MAX_FRAGMENT_LENGTH)
                {
                    hasUpperPosRead = true;
                    strandCount += record.getReadNegativeStrandFlag() ? -1 : 1;
                }

                if(positionsOverlap(varPosMin, varPosMax, readStart, readEnd))
                {
                    SV_LOGGER.trace("var({}) pos({}-{}) read({}-{}) id({}) has ref support",
                            variantPosition, varPosMin, varPosMax, readStart, readEnd, record.getReadName());
                    readSupportsRef = true;
                }
            }

            if(matchesJunction)
                continue;

            if(readSupportsRef)
            {
                ++counts.RefSupport;
            }
            else if(hasLowerPosRead && hasUpperPosRead && strandCount == 0)
            {
                ++counts.RefPairSupport;

                SV_LOGGER.trace("var({}) pos({}-{}) fragment(id={} {}-{}) has ref-pair support",
                        variantPosition, varPosMin, varPosMax, readGroup.id(), readGroupPosMin, readGroupPosMax);
            }

            if(refFragsCap > 0 && counts.total() > refFragsCap)
            {
                SV_LOGGER.debug("var({}:{}) varFrags({}) ref limit({}) reached with support(ref={} pair={})",
                        mChromosome, variantPosition, varFrags, refFragsCap, counts.RefSupport, counts.RefPairSupport);
                break;
            }
        }

        return counts;
    }

    private class SamReadGroup
    {
        public final List<SAMRecord> Reads;

        public SamReadGroup(final SAMRecord record)
        {
            Reads = Lists.newArrayList(record);
        }

        public String id() { return Reads.get(0).getReadName(); }
        public String toString() { return format("id(%s) reads(%d)", id(), Reads.size()); }
    }
}
