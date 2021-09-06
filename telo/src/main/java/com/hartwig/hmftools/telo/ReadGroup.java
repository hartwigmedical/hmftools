package com.hartwig.hmftools.telo;

import static com.hartwig.hmftools.telo.TeloConfig.TE_LOGGER;

import java.util.ArrayList;
import java.util.List;

import com.hartwig.hmftools.common.utils.sv.ChrBaseRegion;

import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMTag;

public class ReadGroup
{
    public enum FragmentType
    {
        UNKNOWN,
        NOT_TELOMERE,
        F1,
        F2,
        F4
    }

    public final List<SAMRecord> Reads = new ArrayList<>();
    public final List<SAMRecord> SupplementaryReads = new ArrayList<>();
    private final String mName;

    public ReadGroup(String name)
    {
        mName = name;
    }

    public final String getName() { return mName; }

    public boolean isComplete()
    {
        if (Reads.isEmpty())
            return false;

        // we check several things
        if (Reads.get(0).getReadPairedFlag() && Reads.size() < 2)
        {
            // we havent got all the reads yet
            return false;
        }

        return findMissingReadBaseRegions().isEmpty();
    }

    public static boolean hasSuppAlignment(final List<SAMRecord> reads)
    {
        return reads.stream().anyMatch(SAMRecord::getSupplementaryAlignmentFlag);
    }

    public boolean contains(SAMRecord read)
    {
        if (!read.getReadName().equals(getName()))
        {
            return false;
        }

        List<SAMRecord> listToLook;
        if (read.getSupplementaryAlignmentFlag())
        {
            listToLook = SupplementaryReads;
        }
        else
        {
            listToLook = Reads;
        }

        // we only check chromosome and alignment start
        return listToLook.stream().anyMatch(x ->
                        x.getFirstOfPairFlag() == read.getFirstOfPairFlag() &&
                        x.getAlignmentStart() == read.getAlignmentStart() &&
                        x.getReferenceName().equals(read.getReferenceName()));
    }

    //public boolean isDuplicate() { return Reads.stream().anyMatch(x -> x.isDuplicate()); }

    public String toString() { return String.format("%s reads(%d) complete(%s)", getName(), Reads.size(), isComplete()); }

    /*
    public String findOtherChromosome(final String chromosome)
    {
        for(ReadRecord read : Reads)
        {
            if(!read.mateChromosome().equals(chromosome))
                return read.mateChromosome();

            if(read.hasSuppAlignment())
                return suppAlignmentChromosome(read.getSuppAlignment());
        }

        return null;
    }*/

    public static final String SUPP_ALIGNMENT_DELIM = ",";

    public static List<ChrBaseRegion> suppAlignmentPositions(final String suppAlignment)
    {
        if(suppAlignment == null)
            return null;

        List<ChrBaseRegion> suppAignPos = new ArrayList<>();

        final String[] supplementaryItems = suppAlignment.split(";");

        for (String si : supplementaryItems)
        {
            final String[] items = si.split(SUPP_ALIGNMENT_DELIM);

            if (items.length < 5)
                return null;

            // supplementary(SA) string attribute looks like
            // 4,191039958,+,68S33M,0,0
            // the first word is the chromosome, the second is the alignment start
            String supplementaryAlignRef = items[0];
            int supplementaryAlignStart = Integer.parseInt(items[1]);
            suppAignPos.add(new ChrBaseRegion(supplementaryAlignRef, supplementaryAlignStart, supplementaryAlignStart));
        }
        return suppAignPos;
    }

    public boolean invariant()
    {
        // check to make sure several things:
        // 1. same record cannot appear more than once
        // 2. Reads cannot contain supplementary
        // 3. SupplementaryReads must only contain supplementary

        if (Reads.size() > 2)
        {
            return false;
        }
        if (!Reads.stream().allMatch(x -> x.getReadName().equals(getName())))
        {
            return false;
        }
        if (!SupplementaryReads.stream().allMatch(x -> x.getReadName().equals(getName())))
        {
            return false;
        }
        if (Reads.stream().anyMatch(SAMRecord::getSupplementaryAlignmentFlag))
        {
            return false;
        }
        if (SupplementaryReads.stream().allMatch(SAMRecord::getSupplementaryAlignmentFlag))
        {
            return false;
        }

        // check if any read has appeared more than once

        return true;
    }

    List<ChrBaseRegion> findMissingReadBaseRegions()
    {
        List<ChrBaseRegion> baseRegions = new ArrayList<>();
        assert(invariant());

        for(SAMRecord read : Reads)
        {
            if(read.getReadPairedFlag() && read.getMateReferenceIndex() != SAMRecord.NO_ALIGNMENT_REFERENCE_INDEX)
            {
                // note that even if the mate unmapped flag is set, we still might not be there
                // unmapped read of a mapped read pair would show up with the mate position
                if(read.getMateReferenceName().isEmpty() || read.getMateAlignmentStart() <= 0)
                {
                    // this shouldn't happen
                    TE_LOGGER.warn("read({}) invalid mate reference, mate ref index({})",
                            read, read.getMateReferenceIndex());
                }
                else if (Reads.stream().noneMatch(
                        x -> x.getFirstOfPairFlag() != read.getFirstOfPairFlag() && // we must check this condition since unmapped read would have
                                // alignment of mate
                                x.getAlignmentStart() == read.getMateAlignmentStart() &&
                                x.getReferenceName().equals(read.getMateReferenceName())))
                {
                    baseRegions.add(new ChrBaseRegion(read.getMateReferenceName(), read.getMateAlignmentStart(), read.getMateAlignmentStart()));
                }
            }

            String saAttribute = read.getStringAttribute(SAMTag.SA.name());
            if (saAttribute != null)
            {
                List<ChrBaseRegion> supplementaryRegions = ReadGroup.suppAlignmentPositions(saAttribute);
                for (ChrBaseRegion br : supplementaryRegions)
                {
                    if (br != null)
                    {
                        // check if this supplementary read exists
                        if (SupplementaryReads.stream()
                                .noneMatch(x -> x.getFirstOfPairFlag() == read.getFirstOfPairFlag() &&
                                        x.getAlignmentStart() == br.start() &&
                                        x.getReferenceName().equals(br.Chromosome)))
                        {
                            baseRegions.add(br);
                        }
                    }
                }
            }
        }
        return baseRegions;
    }
}

