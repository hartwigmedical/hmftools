package com.hartwig.hmftools.esvee.read;

import static com.hartwig.hmftools.esvee.SvConfig.SV_LOGGER;

import java.util.ArrayList;
import java.util.List;

import com.hartwig.hmftools.common.genome.refgenome.RefGenomeInterface;
import com.hartwig.hmftools.esvee.SvConstants;

import org.jetbrains.annotations.Nullable;

import htsjdk.samtools.CigarElement;
import htsjdk.samtools.CigarOperator;

public class RecordNormaliser
{
    private final ReadRescue mReadRescue;
    private final PolyGTrimmer mPolyGTrimmer;

    private final int mSmallIndelMaxEdgeDistance;
    private final int mSmallIndelMinSizeToSoftClip;

    public RecordNormaliser(final RefGenomeInterface referenceGenome)
    {
        mReadRescue = new ReadRescue(referenceGenome);
        mPolyGTrimmer = new PolyGTrimmer(SvConstants.NORMALISERPOLYGLENGTH);

        mSmallIndelMaxEdgeDistance = SvConstants.NORMALISERINDELMAXEDGEDISTANCE;
        mSmallIndelMinSizeToSoftClip = SvConstants.NORMALISERINDELMINSIZETOSOFTCLIP;
    }

    @Nullable
    private CigarElement checkSmallEdgeIndelsForSoftClip(final CigarElement edge, final CigarElement inside, final CigarElement next)
    {
        if(edge.getOperator() != CigarOperator.M || edge.getLength() >= mSmallIndelMaxEdgeDistance)
            return null;
        if(inside.getOperator() != CigarOperator.I && inside.getOperator() != CigarOperator.D)
            return null;
        if(next.getOperator() != CigarOperator.M)
            return null; // Some malformed cigars we see are like: "15M10I100S"
        if(inside.getLength() < mSmallIndelMinSizeToSoftClip)
            return null;
        final int length = inside.getOperator() != CigarOperator.D
                ? edge.getLength() + inside.getLength()
                : edge.getLength();

        return new CigarElement(length, CigarOperator.S);
    }

    /*
    private MutableRecord softclipEdgeIndels(final MutableRecord record)
    {
        if(record.getCigar().getCigarElements().size() < 3)
            return record;

        MutableRecord normalised = record;

        @Nullable
        final CigarElement newStart = checkSmallEdgeIndelsForSoftClip(
                normalised.getCigar().getCigarElement(0),
                normalised.getCigar().getCigarElement(1),
                normalised.getCigar().getCigarElement(2));
        if(newStart != null)
        {
            normalised = record.copyRecord();
            final List<CigarElement> oldElements = normalised.getCigar().getCigarElements();
            final List<CigarElement> newElements = new ArrayList<>();
            newElements.add(newStart);
            for(int i = 2; i < oldElements.size(); i++)
                newElements.add(oldElements.get(i));
            normalised.setAlignmentStart(normalised.getAlignmentStart() + oldElements.get(0).getLength());
            if(oldElements.get(1).getOperator() == CigarOperator.D)
                normalised.setAlignmentStart(normalised.getAlignmentStart() + oldElements.get(1).getLength());
            normalised.setCigar(new Cigar(newElements));
        }

        if(normalised.getCigar().getCigarElements().size() < 3)
            return normalised;

        @Nullable
        final CigarElement newEnd = checkSmallEdgeIndelsForSoftClip(
                normalised.getCigar().getCigarElement(normalised.getCigar().numCigarElements() - 1),
                normalised.getCigar().getCigarElement(normalised.getCigar().numCigarElements() - 2),
                normalised.getCigar().getCigarElement(normalised.getCigar().numCigarElements() - 3));
        if(newEnd != null)
        {
            normalised = record.copyRecord();
            final List<CigarElement> oldElements = normalised.getCigar().getCigarElements();
            final List<CigarElement> newElements = new ArrayList<>();
            for(int i = 0; i < oldElements.size() - 2; i++)
                newElements.add(oldElements.get(i));
            newElements.add(newEnd);
            normalised.setCigar(new Cigar(newElements));
        }

        return normalised;
    }

    private static MutableRecord cleanCigar(final MutableRecord record)
    {
        if(record.getCigar().getCigarElements().size() < 2)
            return record;

        boolean modified = false;
        List<CigarElement> elements = record.getCigar().getCigarElements();
        for(int i = 1; i < elements.size(); i++)
        {
            final CigarElement previous = elements.get(i - 1);
            final CigarElement current = elements.get(i);

            @Nullable
            CigarElement replacement = null;
            if(previous.getOperator() == CigarOperator.S)
            {
                // why would there be 2 soft-clips in a row?
                if(current.getOperator() == CigarOperator.S || current.getOperator() == CigarOperator.I)
                    replacement = new CigarElement(previous.getLength() + current.getLength(), CigarOperator.S);
            }
            if(current.getOperator() == CigarOperator.S)
            {
                if(previous.getOperator() == CigarOperator.S || previous.getOperator() == CigarOperator.I)
                    replacement = new CigarElement(previous.getLength() + current.getLength(), CigarOperator.S);
            }
            if(replacement == null)
                continue;

            if(!modified)
            {
                elements = new ArrayList<>(elements);
                modified = true;
            }

            elements.set(i - 1, replacement);
            elements.remove(i);
            i--;
        }

        if(!modified)
            return record;

        SV_LOGGER.debug("clear cigar and copying read({})", record.getName());

        final MutableRecord normalised = record.copyRecord();
        normalised.setCigar(new Cigar(elements));
        return normalised;
    }
    */

    public void normalise(final Read read)
    {
        // All normalisation atm deals with alignment, so ignore records that are unmapped.
        if(read.isUnmapped())
            return;

        /* CHECK what is actually required and why
        normalised = softclipEdgeIndels(normalised);
        normalised = cleanCigar(normalised);
        normalised = mPolyGTrimmer.trimPolyG(normalised);
        normalised = mReadRescue.rescueRead(normalised);
        */
    }
}
