package com.hartwig.hmftools.redux.common;

import java.util.List;

import com.hartwig.hmftools.common.region.ChrBaseRegion;
import com.hartwig.hmftools.common.bam.SupplementaryReadData;

import htsjdk.samtools.SAMRecord;

public enum FilterReadsType
{
    NONE,
    READ, // ignore a read if its coords are outside the specified region
    MATE, // ignore a read if its mate is outside the specified regions
    MATE_AND_SUPP; // as above but also check supplementaries


    public boolean filterMates() { return this == MATE || this == MATE_AND_SUPP; }
    public boolean filterSupplementaries() { return this == MATE_AND_SUPP; }

    public static boolean readOutsideSpecifiedRegions(
            final SAMRecord read, final List<ChrBaseRegion> regions, final List<String> chromosomes, final FilterReadsType filterType)
    {
        if(filterType == NONE)
            return false;

        if(!chromosomes.isEmpty())
        {
            // ignore the chromosome since reads are sliced by chromosome, and this otherwise excludes unmapped reads which
            // do not have their contig set (although reference name is set)

            //if(chromosomes.stream().noneMatch(x -> x.equals(read.getReferenceName())))
            //    return true;

            // any mates or supplementaries must also be within the regions specified
            if(filterType.filterMates() && chromosomes.stream().noneMatch(x -> x.equals(read.getMateReferenceName())))
                return true;
        }

        if(!regions.isEmpty())
        {
            if(filterType == READ && regions.stream().noneMatch(x -> x.containsPosition(read.getReferenceName(), read.getAlignmentStart())))
                return true;

            // any mates or supplementaries must also be within the regions specified
            if(filterType.filterMates() && regions.stream().noneMatch(x -> x.containsPosition(read.getMateReferenceName(), read.getMateAlignmentStart())))
                return true;
        }

        // by default ignore checking supplementaries since a) they aren't marked as duplicates by other tools and b) they shouldn't be
        // a reason to ignore a primary read since that then impacts duplicate classification
        if(filterType.filterSupplementaries())
        {
            SupplementaryReadData suppData = SupplementaryReadData.extractAlignment(read);

            if(suppData != null)
            {
                if(!regions.isEmpty() && regions.stream().noneMatch(x -> x.containsPosition(suppData.Chromosome, suppData.Position)))
                    return true;

                if(!chromosomes.isEmpty() && chromosomes.stream().noneMatch(x -> x.equals(suppData.Chromosome)))
                    return true;
            }
        }

        return false;
    }

}
