package com.hartwig.hmftools.svprep.depth;

import static java.lang.String.format;

import static com.hartwig.hmftools.common.sv.StructuralVariantFactory.REFERENCE_BREAKEND_READPAIR_COVERAGE;
import static com.hartwig.hmftools.common.sv.StructuralVariantFactory.REFERENCE_BREAKEND_READ_COVERAGE;
import static com.hartwig.hmftools.common.utils.sv.BaseRegion.positionWithin;
import static com.hartwig.hmftools.common.utils.sv.SvCommonUtils.NEG_ORIENT;
import static com.hartwig.hmftools.common.utils.sv.SvCommonUtils.POS_ORIENT;
import static com.hartwig.hmftools.svprep.SvCommon.SV_LOGGER;
import static com.hartwig.hmftools.svprep.SvConstants.MAX_FRAGMENT_LENGTH;

import java.io.File;
import java.util.List;
import java.util.Map;
import java.util.concurrent.Callable;

import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import com.hartwig.hmftools.common.samtools.BamSlicer;
import com.hartwig.hmftools.common.utils.sv.ChrBaseRegion;

import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;
import htsjdk.variant.variantcontext.VariantContext;

public class DepthTask implements Callable
{
    private final List<VariantContext> mVariantsList;
    private final String mChromosome;
    private int mBamRecords;

    private final SamReader mSamReader;
    private final BamSlicer mBamSlicer;

    private final Map<String,SamReadGroup> mReadGroups;

    public DepthTask(final String chromosome, final String refGenomeFile, final String bamFile)
    {
        mVariantsList = Lists.newArrayList();
        mChromosome = chromosome;
        mBamRecords = 0;

        mSamReader = SamReaderFactory.makeDefault().referenceSequence(new File(refGenomeFile)).open(new File(bamFile));
        mBamSlicer = new BamSlicer(0, true, true, true);
        mReadGroups = Maps.newHashMap();
    }

    public List<VariantContext> variants() { return mVariantsList; }

    @Override
    public Long call()
    {
        SV_LOGGER.info("chr({}) processing {} variants", mChromosome, mVariantsList.size());

        int processed = 0;
        for(VariantContext variant : mVariantsList)
        {
            retrieveDepth(variant);
            ++processed;

            if((processed % 100) == 0)
            {
                SV_LOGGER.debug("chr({}) processed {} variants", mChromosome, processed);
            }
        }

        SV_LOGGER.info("chr({}) complete for {} variants, total reads({})", mChromosome, processed, mBamRecords);

        return (long)0;
    }

    private void retrieveDepth(final VariantContext variant)
    {
        // retrieve the depth to set these 2 values in the VCF:
        // REF = reads aligned directly over the breakend (excluding reads that would support the junction)
        // REFPAIR = fragments with 1 read aligned to the left and 1 read aligned to the right of the breakend with proper orientaition and
        // FragmentSize < max size from distribution
        mReadGroups.clear();

        int variantPosition = variant.getStart();

        ChrBaseRegion region = new ChrBaseRegion(
                variant.getContig(), variantPosition - MAX_FRAGMENT_LENGTH, variantPosition + MAX_FRAGMENT_LENGTH);

        mBamSlicer.slice(mSamReader, Lists.newArrayList(region), this::processSamRecord);

        int refSupport = 0;
        int refPairSupport = 0;

        for(SamReadGroup readGroup : mReadGroups.values())
        {
            boolean readSupportsRef = false;
            boolean hasLowerPosRead = false;
            boolean hasUpperPosRead = false;

            for(SAMRecord record : readGroup.Reads)
            {
                int readStart = record.getAlignmentStart();
                int readEnd = record.getAlignmentEnd();

                byte orientation = !record.getReadNegativeStrandFlag() ? POS_ORIENT : NEG_ORIENT;

                if(orientation == POS_ORIENT && readEnd < variantPosition)
                {
                    hasLowerPosRead = true;
                }
                else if(orientation == NEG_ORIENT && readStart > variantPosition)
                {
                    hasUpperPosRead = true;
                }

                if(positionWithin(variantPosition, readStart + 1, readEnd - 1))
                {
                    readSupportsRef = true;
                }
            }

            if(readSupportsRef)
                ++refSupport;
            else if(hasLowerPosRead && hasUpperPosRead)
                ++refPairSupport;
        }

        int vcfRef = variant.getAttributeAsInt(REFERENCE_BREAKEND_READ_COVERAGE, 0);
        int vcfRefPair = variant.getAttributeAsInt(REFERENCE_BREAKEND_READPAIR_COVERAGE, 0);

        SV_LOGGER.debug("var({}:{}) vcf(REF={} REFPAIR={}) calc(REF={} REFPAIR={}",
                variant.getContig(), variantPosition, vcfRef, vcfRefPair, refSupport, refPairSupport);
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
