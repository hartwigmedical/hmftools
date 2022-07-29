package com.hartwig.hmftools.svprep.depth;

import static java.lang.Math.abs;
import static java.lang.Math.max;
import static java.lang.String.format;

import static com.hartwig.hmftools.common.sv.ExcludedRegions.getPolyGRegions;
import static com.hartwig.hmftools.common.sv.StructuralVariantFactory.IHOMPOS;
import static com.hartwig.hmftools.common.sv.StructuralVariantFactory.REFERENCE_BREAKEND_READPAIR_COVERAGE;
import static com.hartwig.hmftools.common.sv.StructuralVariantFactory.REFERENCE_BREAKEND_READ_COVERAGE;
import static com.hartwig.hmftools.common.sv.StructuralVariantFactory.isSingleBreakend;
import static com.hartwig.hmftools.common.utils.sv.BaseRegion.positionWithin;
import static com.hartwig.hmftools.common.utils.sv.BaseRegion.positionsOverlap;
import static com.hartwig.hmftools.common.utils.sv.SvCommonUtils.NEG_ORIENT;
import static com.hartwig.hmftools.common.utils.sv.SvCommonUtils.POS_ORIENT;
import static com.hartwig.hmftools.svprep.SvCommon.SV_LOGGER;
import static com.hartwig.hmftools.svprep.SvConstants.MAX_FRAGMENT_LENGTH;
import static com.hartwig.hmftools.svprep.depth.DepthAnnotator.VCF_TAG_REFPAIR_GRIDSS;
import static com.hartwig.hmftools.svprep.depth.DepthAnnotator.VCF_TAG_REF_GRIDSS;

import java.io.File;
import java.util.List;
import java.util.Map;
import java.util.concurrent.Callable;
import java.util.stream.Collectors;

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
import htsjdk.variant.variantcontext.VariantContextBuilder;

public class DepthTask implements Callable
{
    private final DepthConfig mConfig;
    private final Map<String,Integer> mSampleVcfGenotypeIds;
    private final List<VariantContext> mVariantsList;
    private final List<VariantContext> mNewVariantsList;
    private final String mChromosome;
    private int mBamRecords;

    private final List<SamReader> mSamReaders;
    private final BamSlicer mBamSlicer;

    private final Map<String,SamReadGroup> mReadGroups;

    private final List<ChrBaseRegion> mExcludedRegions;
    private final PerformanceCounter mPerfCounter;

    public DepthTask(final String chromosome, final DepthConfig config, final Map<String,Integer> sampleVcfGenotypeIds)
    {
        mConfig = config;
        mChromosome = chromosome;
        mSampleVcfGenotypeIds = sampleVcfGenotypeIds;

        mVariantsList = Lists.newArrayList();
        mNewVariantsList = Lists.newArrayList();
        mBamRecords = 0;

        mBamSlicer = new BamSlicer(0, false, true, false);
        mReadGroups = Maps.newHashMap();

        mSamReaders = Lists.newArrayList();

        for(String bamFile : mConfig.BamFiles)
        {
            mSamReaders.add(SamReaderFactory.makeDefault().referenceSequence(new File(mConfig.RefGenome)).open(new File(bamFile)));
        }

        mExcludedRegions = getPolyGRegions(mConfig.RefGenVersion).stream()
                .filter(x -> x.Chromosome.equals(chromosome)).collect(Collectors.toList());

        mPerfCounter = new PerformanceCounter("Slice");
    }

    public String chromosome() { return mChromosome; }
    public void addVariants(final List<VariantContext> variants) { mVariantsList.addAll(variants); }
    public List<VariantContext> newVariants() { return mNewVariantsList; }
    public PerformanceCounter getPerfCounter() { return mPerfCounter; }

    @Override
    public Long call()
    {
        SV_LOGGER.info("chr({}) processing {} variants", mChromosome, mVariantsList.size());

        int processed = 0;
        for(VariantContext variant : mVariantsList)
        {
            retrieveDepth(variant);
            ++processed;

            if((processed % 1000) == 0)
            {
                SV_LOGGER.debug("chr({}) processed {} variants", mChromosome, processed);
            }
        }

        SV_LOGGER.info("chr({}) complete for {} variants, total reads({})", mChromosome, processed, mBamRecords);
        System.gc();

        return (long)0;
    }

    private void retrieveDepth(final VariantContext variant)
    {
        mPerfCounter.start();

        // retrieve the depth to set these 2 values in the VCF:
        // REF = reads aligned directly over the breakend (excluding reads that would support the junction)
        // REFPAIR = fragments with 1 read aligned to the left and 1 read aligned to the right of the breakend with proper orientation and
        // FragmentSize < max size from distribution
        int variantPosition = variant.getStart();
        byte variantOrientation = getOrientation(variant);

        ChrBaseRegion region = new ChrBaseRegion(
                variant.getContig(), variantPosition - MAX_FRAGMENT_LENGTH, variantPosition + MAX_FRAGMENT_LENGTH);

        final int[] homology = {0, 0};

        if(variant.hasAttribute(IHOMPOS))
        {
            final List<Integer> ihompos = variant.getAttributeAsIntList(IHOMPOS, 0);
            homology[0] = ihompos.get(0);
            homology[1] = ihompos.get(1);
        }

        int varPosStart = variantPosition + homology[0];
        int varPosEnd = variantPosition + homology[1];

        VariantContext newVariant = new VariantContextBuilder(variant)
                .genotypes(variant.getGenotypes())
                .filters(variant.getFilters())
                .make();

        RefSupportCounts totalCounts = new RefSupportCounts();

        for(int i = 0; i < mConfig.Samples.size(); ++i)
        {
            String sampleId = mConfig.Samples.get(i);
            SamReader samReader = mSamReaders.get(i);

            mReadGroups.clear();

            if(mExcludedRegions.stream().anyMatch(x -> x.containsPosition(variantPosition)))
            {
                // decide how to populate
                SV_LOGGER.info("sample({}) var({}:{}) skipped in excluded region", sampleId, variant.getContig(), variantPosition);
                continue;
            }

            mBamSlicer.slice(samReader, Lists.newArrayList(region), this::processSamRecord);

            RefSupportCounts sampleCounts = calculateSupport(variantPosition, varPosStart, varPosEnd, variantOrientation);
            totalCounts.RefSupport += sampleCounts.RefSupport;
            totalCounts.RefPairSupport += sampleCounts.RefPairSupport;

            int vcfSampleIndex = mSampleVcfGenotypeIds.get(sampleId);
            Genotype genotype = newVariant.getGenotype(vcfSampleIndex);

            setRefDepthValue(genotype, sampleCounts.RefSupport, REFERENCE_BREAKEND_READ_COVERAGE, VCF_TAG_REF_GRIDSS);
            setRefDepthValue(genotype, sampleCounts.RefPairSupport, REFERENCE_BREAKEND_READPAIR_COVERAGE, VCF_TAG_REFPAIR_GRIDSS);

            int vcfRef = getGenotypeAttributeAsInt(genotype, REFERENCE_BREAKEND_READ_COVERAGE, 0);
            int vcfRefPair = getGenotypeAttributeAsInt(genotype, REFERENCE_BREAKEND_READPAIR_COVERAGE, 0);

            boolean hasDiffs = hasDiff(vcfRef, sampleCounts.RefSupport) || hasDiff(vcfRefPair, sampleCounts.RefPairSupport);

            Level level = hasDiffs ? Level.DEBUG : Level.TRACE;

            SV_LOGGER.log(level, "sample({}) var({}:{}) {} REF(vcf={} calc={}) REFPAIR(vcf={} calc={})",
                    sampleId, variant.getContig(), variantPosition, hasDiffs ? "has diffs" : "equal",
                    vcfRef, sampleCounts.RefSupport, vcfRefPair, sampleCounts.RefPairSupport);
        }

        setRefDepthValue(newVariant, totalCounts.RefSupport, REFERENCE_BREAKEND_READ_COVERAGE, VCF_TAG_REF_GRIDSS);
        setRefDepthValue(newVariant, totalCounts.RefPairSupport, REFERENCE_BREAKEND_READPAIR_COVERAGE, VCF_TAG_REFPAIR_GRIDSS);

        mNewVariantsList.add(newVariant);

        mPerfCounter.stop();

        if(mPerfCounter.getLastTime() > 2)
        {
            SV_LOGGER.warn("var({}:{}) high depth retrieval time({}",
                    variant.getContig(), variantPosition, format("%.3f", mPerfCounter.getLastTime()));
        }
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

    private static final int ABS_DIFF_MAX = 5;
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
    }

    private RefSupportCounts calculateSupport(int variantPosition, int varPosStart, int varPosEnd, byte variantOrientation)
    {
        RefSupportCounts counts = new RefSupportCounts();

        for(SamReadGroup readGroup : mReadGroups.values())
        {
            boolean readSupportsRef = false;
            boolean hasLowerPosRead = false;
            boolean hasUpperPosRead = false;
            boolean matchesJunction = false;

            boolean isValidSpanningPair = readGroup.Reads.size() == 2
                    && abs(readGroup.Reads.get(0).getInferredInsertSize()) < MAX_FRAGMENT_LENGTH
                    && readGroup.Reads.get(0).getReadNegativeStrandFlag() != readGroup.Reads.get(1).getReadNegativeStrandFlag();

            for(SAMRecord record : readGroup.Reads)
            {
                int readStart = record.getAlignmentStart();
                int readEnd = record.getAlignmentEnd();

                if((variantOrientation == NEG_ORIENT && positionWithin(readStart, varPosStart, varPosEnd))
                || (variantOrientation == POS_ORIENT && positionWithin(readEnd, varPosStart, varPosEnd))) // exact SC match
                {
                    SV_LOGGER.trace("var({}) pos({}-{}) read({}-{}) id({}) at junction",
                            variantPosition, varPosStart, varPosEnd, readStart, readEnd, record.getReadName());
                    matchesJunction = true;
                    break;
                }

                byte orientation = !record.getReadNegativeStrandFlag() ? POS_ORIENT : NEG_ORIENT;

                if(orientation == POS_ORIENT && readEnd < varPosStart)
                {
                    hasLowerPosRead = true;
                }
                else if(orientation == NEG_ORIENT && readStart > varPosEnd)
                {
                    hasUpperPosRead = true;
                }

                if(positionsOverlap(varPosStart, varPosEnd, readStart, readEnd))
                {
                    SV_LOGGER.trace("var({}) pos({}-{}) read({}-{}) id({}) supports",
                            variantPosition, varPosStart, varPosEnd, readStart, readEnd, record.getReadName());
                    readSupportsRef = true;
                }
            }

            if(!matchesJunction)
            {
                if(readSupportsRef)
                    ++counts.RefSupport;
                else if(isValidSpanningPair && hasLowerPosRead && hasUpperPosRead)
                    ++counts.RefPairSupport;
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
