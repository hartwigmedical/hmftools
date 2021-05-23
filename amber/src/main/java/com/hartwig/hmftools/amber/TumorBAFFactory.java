package com.hartwig.hmftools.amber;

import static com.hartwig.hmftools.common.amber.BaseDepthFactory.getBaseQuality;
import static com.hartwig.hmftools.common.amber.BaseDepthFactory.indel;

import com.hartwig.hmftools.common.amber.BaseDepth;
import com.hartwig.hmftools.common.amber.ModifiableTumorBAF;

import org.jetbrains.annotations.NotNull;

import htsjdk.samtools.SAMRecord;

class TumorBAFFactory
{
    private final int mMinBaseQuality;

    TumorBAFFactory(final int minBaseQuality)
    {
        mMinBaseQuality = minBaseQuality;
    }

    @NotNull
    public static ModifiableTumorBAF create(@NotNull final BaseDepth normal)
    {
        return ModifiableTumorBAF.create()
                .from(normal)
                .setRef(normal.ref().toString())
                .setAlt(normal.alt().toString())
                .setNormalReadDepth(normal.readDepth())
                .setNormalRefSupport(normal.refSupport())
                .setNormalAltSupport(normal.altSupport())
                .setTumorIndelCount(0)
                .setTumorReadDepth(0)
                .setTumorRefSupport(0)
                .setTumorAltQuality(0)
                .setTumorAltSupport(0);
    }

    void addEvidence(@NotNull final ModifiableTumorBAF evidence, @NotNull final SAMRecord samRecord)
    {
        int quality = getBaseQuality(evidence, samRecord);
        if(quality >= mMinBaseQuality)
        {
            evidence.setTumorReadDepth(evidence.tumorReadDepth() + 1);
            int bafPosition = (int) evidence.position();
            int readPosition = samRecord.getReadPositionAtReferencePosition(bafPosition);
            if(readPosition != 0)
            {
                if(!indel(bafPosition, readPosition, samRecord))
                {
                    final String base = String.valueOf(samRecord.getReadString().charAt(readPosition - 1));
                    if(base.equals(evidence.ref()))
                    {
                        evidence.setTumorRefSupport(evidence.tumorRefSupport() + 1);
                    }
                    else if(base.equals(evidence.alt()))
                    {
                        evidence.setTumorAltSupport(evidence.tumorAltSupport() + 1);
                        evidence.setTumorAltQuality(evidence.tumorAltQuality() + quality);
                    }
                }
                else
                {
                    evidence.setTumorIndelCount(evidence.tumorIndelCount() + 1);
                }
            }
        }
    }
}
