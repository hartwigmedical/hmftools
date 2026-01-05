package com.hartwig.hmftools.sage.filter;

import static java.lang.Math.abs;

import static com.hartwig.hmftools.common.bam.SamRecordUtils.mateNegativeStrand;
import static com.hartwig.hmftools.common.region.BaseRegion.positionWithin;
import static com.hartwig.hmftools.common.region.BaseRegion.positionsOverlap;
import static com.hartwig.hmftools.sage.SageConstants.CHIMERIC_FRAGMENT_LENGTH_MAX;
import static com.hartwig.hmftools.sage.SageConstants.NEAR_INDEL_MIN_VAF;
import static com.hartwig.hmftools.sage.SageConstants.NEAR_INDEL_PROXIMITY;
import static com.hartwig.hmftools.sage.filter.SoftFilter.TUMOR_FILTERS;

import java.util.List;

import com.google.common.annotations.VisibleForTesting;
import com.hartwig.hmftools.sage.common.SageVariant;

import htsjdk.samtools.SAMRecord;

public final class FilterUtils
{
    public static boolean isChimericRead(final SAMRecord record)
    {
        if(record.getReadPairedFlag())
        {
            if(record.getMateUnmappedFlag())
                return true;

            // inter-chromosomal
            if(!record.getReferenceName().equals(record.getMateReferenceName()))
                return true;

            // inversion
            if(record.getReadNegativeStrandFlag() == mateNegativeStrand(record))
                return true;
        }

        // or a fragment length outside the expected maximum
        if(abs(record.getInferredInsertSize()) > CHIMERIC_FRAGMENT_LENGTH_MAX)
            return true;

        return false;
    }

    public static void setNearByIndelStatusPreFilter(final List<SageVariant> sageVariants)
    {
        // look forward and backwards from this indel and mark other variants which fall within its bounds
        for(int index = 0; index < sageVariants.size(); ++index)
        {
            SageVariant variant = sageVariants.get(index);

            if(!variant.isIndel() || variant.variant().indelLengthAbs() < 2)
                continue;

            double af = variant.tumorReadCounters().get(0).vaf();

            if(af < NEAR_INDEL_MIN_VAF)
                continue;

            int indelStart = variant.position() - NEAR_INDEL_PROXIMITY;
            int indelEnd = variant.variant().positionEnd() + NEAR_INDEL_PROXIMITY;

            for(int i = 0; i <= 1; ++i)
            {
                boolean searchUp = (i == 0);

                int otherIndex = searchUp ? index + 1 : index - 1;

                while(otherIndex >= 0 && otherIndex < sageVariants.size())
                {
                    SageVariant otherVar = sageVariants.get(otherIndex);

                    int otherPosEnd = otherVar.variant().positionEnd();

                    if(positionsOverlap(indelStart, indelEnd, otherVar.position(), otherPosEnd))
                    {
                        otherVar.setNearMultiBaseIndel();
                    }
                    else
                    {
                        if(searchUp && otherVar.position() > indelEnd)
                            break;
                        else if(!searchUp && otherPosEnd < indelStart)
                            break;
                    }

                    if(searchUp)
                        ++otherIndex;
                    else
                        --otherIndex;
                }
            }
        }
    }

    @VisibleForTesting
    public static void setNearByIndelStatusPostFilter(final List<SageVariant> sageVariants)
    {
        // look forward and backwards from this indel and mark other variants which fall within its bounds
        for(int index = 0; index < sageVariants.size(); ++index)
        {
            SageVariant variant = sageVariants.get(index);

            if(!variant.isIndel())
                continue;

            // ignore if filtered other than by germline-only filters
            if(!variant.isPassing() && variant.filters().stream().anyMatch(TUMOR_FILTERS::contains))
                continue;

            for(int i = 0; i <= 1; ++i)
            {
                boolean searchUp = (i == 0);

                int otherIndex = searchUp ? index + 1 : index - 1;

                while(otherIndex >= 0 && otherIndex < sageVariants.size())
                {
                    SageVariant otherVar = sageVariants.get(otherIndex);

                    if(positionWithin(otherVar.position(), variant.readContext().AlignmentStart, variant.readContext().AlignmentEnd))
                    {
                        otherVar.setNearIndel();
                    }
                    else
                    {
                        if(searchUp && otherVar.position() > variant.readContext().AlignmentEnd)
                            break;
                        else if(!searchUp && otherVar.position() < variant.readContext().AlignmentStart)
                            break;
                    }

                    if(searchUp)
                        ++otherIndex;
                    else
                        --otherIndex;
                }
            }
        }
    }
}
