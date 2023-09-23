package com.hartwig.hmftools.common.region;

import static com.hartwig.hmftools.common.region.BaseRegion.checkMergeOverlaps;
import static com.hartwig.hmftools.common.region.BaseRegion.positionWithin;
import static com.hartwig.hmftools.common.region.BaseRegion.positionsOverlap;
import static com.hartwig.hmftools.common.utils.file.CommonFields.FLD_CHROMOSOME;
import static com.hartwig.hmftools.common.utils.file.CommonFields.FLD_POSITION_END;
import static com.hartwig.hmftools.common.utils.file.CommonFields.FLD_POSITION_START;
import static com.hartwig.hmftools.common.utils.file.CommonFields.FLD_POS_END;
import static com.hartwig.hmftools.common.utils.file.CommonFields.FLD_POS_START;
import static com.hartwig.hmftools.common.utils.sv.StartEndIterator.SE_END;
import static com.hartwig.hmftools.common.utils.sv.StartEndIterator.SE_START;

import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.List;
import java.util.Map;
import java.util.NoSuchElementException;

import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import com.hartwig.hmftools.common.genome.chromosome.HumanChromosome;
import com.hartwig.hmftools.common.genome.region.GenomeRegion;
import com.hartwig.hmftools.common.genome.region.GenomeRegions;
import com.hartwig.hmftools.common.utils.file.FileDelimiters;
import com.hartwig.hmftools.common.utils.file.FileReaderUtils;

import org.jetbrains.annotations.NotNull;

public class ChrBaseRegion implements Cloneable, Comparable<ChrBaseRegion>
{
    public final String Chromosome;
    private int mStart;
    private int mEnd;

    public ChrBaseRegion(final String chromosome, final int[] positions)
    {
        Chromosome = chromosome;
        mStart = positions[SE_START];
        mEnd = positions[SE_END];
    }

    public ChrBaseRegion(final String chromosome, final int posStart, final int posEnd)
    {
        Chromosome = chromosome;
        mStart = posStart;
        mEnd = posEnd;
    }

    public static ChrBaseRegion from(final GenomeRegion region) { return new ChrBaseRegion(region.chromosome(), region.start(), region.end()); }
    public GenomeRegion genomeRegion() { return GenomeRegions.create(chromosome(), start(), end()); }

    public int start() { return mStart; }
    public int end() { return mEnd; }

    public int position(int which)
    {
        if(which == SE_START)
        {
            return mStart;
        }
        else if(which == SE_END)
        {
            return mEnd;
        }
        throw new NoSuchElementException();
    }

    public String chromosome() { return Chromosome; }
    
    public void setStart(int pos) { mStart = pos; }
    public void setEnd(int pos) { mEnd = pos; }

    public int baseLength() { return length() + 1; }
    public int length() { return mEnd - mStart; }

    public boolean isValid() { return HumanChromosome.contains(Chromosome) && hasValidPositions(); }
    public boolean hasValidPositions() { return mStart > 0 & mEnd >= mStart; }

    public boolean overlaps(final ChrBaseRegion other)
    {
        if(!Chromosome.equals(other.Chromosome))
            return false;

        return positionsOverlap(mStart, mEnd, other.mStart, other.mEnd);
    }

    public boolean containsPosition(int position) { return positionWithin(position, start(), end()); }

    public boolean containsPosition(final String chromosome, int position)
    {
        return Chromosome.equals(chromosome) && positionWithin(position, start(), end());
    }

    public boolean matches(final ChrBaseRegion other)
    {
        return Chromosome.equals(other.Chromosome) && start() == other.start() && end() == other.end();
    }

    public static boolean containsPosition(final List<ChrBaseRegion> regions, final String chromosome, final int position)
    {
        return regions.stream().anyMatch(x -> x.containsPosition(chromosome, position));
    }

    public static boolean overlaps(final List<ChrBaseRegion> regions, final ChrBaseRegion region)
    {
        return regions.stream().anyMatch(x -> x.overlaps(region));
    }

    public String toString() { return String.format("%s:%d-%d", Chromosome, mStart, mEnd); }

    @Override
    public ChrBaseRegion clone()
    {
        try
        {
            ChrBaseRegion br = (ChrBaseRegion) super.clone();
            br.mStart = mStart;
            br.mEnd = mEnd;
            return br;
        }
        catch (CloneNotSupportedException e)
        {
            // Will not happen in this case
            return null;
        }
    }

    @Override
    public boolean equals(Object obj)
    {
        if(obj == this)
            return true;

        if(obj == null)
            return false;

        if(!getClass().equals(obj.getClass()))
            return false;

        ChrBaseRegion other = (ChrBaseRegion) obj;
        return matches(other);
    }

    @Override
    public int hashCode()
    {
        int result = 31 + Chromosome.hashCode();
        result = 31 * result + mStart;
        result = 31 * result + mEnd;
        return result;
    }

    @Override
    public int compareTo(@NotNull final ChrBaseRegion other)
    {
        if(Chromosome.equals(other.Chromosome))
        {
            if(start() < other.start())
            {
                return -1;
            }
            else if(start() == other.start())
            {
                return 0;
            }
            return 1;
        }

        int rank1 = HumanChromosome.chromosomeRank(Chromosome);
        int rank2 = HumanChromosome.chromosomeRank(other.Chromosome);

        if(rank1 == rank2)
            return 0;

        return rank1 < rank2 ? -1 : 1;
    }

    public static Map<String,List<BaseRegion>> loadChrBaseRegions(final String filename)
    {
        return loadChrBaseRegions(filename, false);
    }

    public static Map<String,List<BaseRegion>> loadChrBaseRegions(final String filename, boolean isBedFile)
    {
        Map<String,List<BaseRegion>> chrRegionsMap = Maps.newHashMap();

        if(filename == null)
            return chrRegionsMap;

        try
        {
            List<String> lines = Files.readAllLines(Paths.get(filename));
            String delim = FileDelimiters.inferFileDelimiter(filename);

            String header = lines.get(0);
            lines.remove(0);
            Map<String,Integer> fieldIndexMap = FileReaderUtils.createFieldsIndexMap(header, delim);

            int chrIndex = fieldIndexMap.get(FLD_CHROMOSOME);

            Integer posStartIndex = fieldIndexMap.containsKey(FLD_POSITION_START) ?
                    fieldIndexMap.get(FLD_POSITION_START) : fieldIndexMap.get(FLD_POS_START);

            Integer posEndIndex = fieldIndexMap.containsKey(FLD_POSITION_END) ?
                    fieldIndexMap.get(FLD_POSITION_END) : fieldIndexMap.get(FLD_POS_END);

            for(String line : lines)
            {
                final String[] values = line.split(delim, -1);

                String chromosome = values[chrIndex];
                int posStart = Integer.parseInt(values[posStartIndex]);

                if(isBedFile)
                    ++posStart;

                int posEnd = Integer.parseInt(values[posEndIndex]);

                List<BaseRegion> regions = chrRegionsMap.get(chromosome);

                if(regions == null)
                {
                    regions = Lists.newArrayList();
                    chrRegionsMap.put(chromosome, regions);
                }

                regions.add(new BaseRegion(posStart, posEnd));
            }

            chrRegionsMap.values().forEach(x -> checkMergeOverlaps(x));

            return chrRegionsMap;
        }
        catch(IOException e)
        {
            return null;
        }
    }
}

