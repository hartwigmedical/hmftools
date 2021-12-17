package com.hartwig.hmftools.common.amber;

import com.hartwig.hmftools.common.genome.position.GenomePosition;
import com.hartwig.hmftools.common.samtools.SamRecordUtils;

import org.jetbrains.annotations.NotNull;

import htsjdk.samtools.SAMRecord;
import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.Genotype;
import htsjdk.variant.variantcontext.VariantContext;

public class BaseDepthFactory
{
    @NotNull
    public static ModifiableBaseDepth create(@NotNull final BaseDepth pos)
    {
        return ModifiableBaseDepth.create().from(pos).setIndelCount(0).setRefSupport(0).setAltSupport(0).setReadDepth(0);
    }

    @NotNull
    public static BaseDepth fromVariantContext(@NotNull final VariantContext context)
    {
        final Allele refAllele = context.getReference();
        final Allele altAllele = context.getAlternateAllele(0);
        final Genotype normal = context.getGenotype(0);

        return ModifiableBaseDepth.create()
                .setChromosome(context.getContig())
                .setPosition(context.getStart())
                .setRef(BaseDepth.Base.valueOf(refAllele.getBaseString()))
                .setAlt(BaseDepth.Base.valueOf(altAllele.getBaseString()))
                .setAltSupport(normal.getAD()[0])
                .setRefSupport(normal.getAD()[1])
                .setIndelCount(0)
                .setReadDepth(normal.getDP());
    }

    @NotNull
    public static ModifiableBaseDepth fromAmberSite(@NotNull final AmberSite site)
    {
        return ModifiableBaseDepth.create()
                .setChromosome(site.chromosome())
                .setPosition(site.position())
                .setRef(BaseDepth.Base.valueOf(site.ref()))
                .setAlt(BaseDepth.Base.valueOf(site.alt()))
                .setAltSupport(0)
                .setRefSupport(0)
                .setIndelCount(0)
                .setReadDepth(0);
    }

    private final int minBaseQuality;

    public BaseDepthFactory(final int minBaseQuality)
    {
        this.minBaseQuality = minBaseQuality;
    }

    public void addEvidence(@NotNull final ModifiableBaseDepth evidence, @NotNull final SAMRecord samRecord)
    {
        int quality = getBaseQuality(evidence, samRecord);
        if(quality >= minBaseQuality)
        {
            evidence.setReadDepth(evidence.readDepth() + 1);

            int bafPosition = evidence.position();
            int readPosition = samRecord.getReadPositionAtReferencePosition(bafPosition);
            if(readPosition != 0)
            {
                if(!indel(bafPosition, readPosition, samRecord))
                {
                    final char baseChar = samRecord.getReadString().charAt(readPosition - 1);
                    final BaseDepth.Base base = BaseDepth.Base.valueOf(String.valueOf(baseChar).toUpperCase());
                    if(base.equals(evidence.ref()))
                    {
                        evidence.setRefSupport(evidence.refSupport() + 1);
                    }
                    else if(base.equals(evidence.alt()))
                    {
                        evidence.setAltSupport(evidence.altSupport() + 1);
                    }
                }
                else
                {
                    evidence.setIndelCount(evidence.indelCount() + 1);
                }
            }
        }
    }

    public static boolean indel(int bafPosition, int readPosition, @NotNull final SAMRecord samRecord)
    {
        if(samRecord.getAlignmentEnd() > bafPosition)
        {
            // Delete?
            if(samRecord.getReadPositionAtReferencePosition(bafPosition + 1) == 0)
            {
                return true;
            }

            // Insert?
            return samRecord.getReferencePositionAtReadPosition(readPosition + 1) != bafPosition + 1;
        }

        return false;
    }

    public static int getBaseQuality(@NotNull final GenomePosition position, @NotNull final SAMRecord samRecord)
    {
        // Get quality of base after del if necessary
        for(int pos = position.position(); pos <= samRecord.getAlignmentEnd(); pos++)
        {
            int readPosition = samRecord.getReadPositionAtReferencePosition(pos);
            if(readPosition != 0)
            {
                return SamRecordUtils.getBaseQuality(samRecord, readPosition);
            }
        }

        return 0;
    }
}
