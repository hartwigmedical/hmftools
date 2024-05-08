package com.hartwig.hmftools.common.region;

import static com.hartwig.hmftools.common.region.BaseRegion.positionWithin;
import static com.hartwig.hmftools.common.region.BaseRegion.positionsOverlap;
import static com.hartwig.hmftools.common.utils.file.CommonFields.FLD_CHROMOSOME;
import static com.hartwig.hmftools.common.utils.file.CommonFields.FLD_POSITION;
import static com.hartwig.hmftools.common.utils.file.CommonFields.FLD_POSITION_END;
import static com.hartwig.hmftools.common.utils.file.CommonFields.FLD_POSITION_START;
import static com.hartwig.hmftools.common.utils.file.CommonFields.FLD_POS_END;
import static com.hartwig.hmftools.common.utils.file.CommonFields.FLD_POS_START;
import static com.hartwig.hmftools.common.utils.file.FileWriterUtils.createBufferedReader;
import static com.hartwig.hmftools.common.utils.sv.StartEndIterator.SE_END;
import static com.hartwig.hmftools.common.utils.sv.StartEndIterator.SE_START;

import java.io.BufferedReader;
import java.io.IOException;
import java.util.Collections;
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

    public BaseRegion baseRegion()
    {
        return new BaseRegion(mStart, mEnd);
    }

    public boolean overlaps(final ChrBaseRegion other)
    {
        if(!Chromosome.equals(other.Chromosome))
            return false;

        return positionsOverlap(mStart, mEnd, other.mStart, other.mEnd);
    }

    public boolean overlaps(final String chromosome, final int posStart, final int posEnd)
    {
        return Chromosome.equals(chromosome) && positionsOverlap(mStart, mEnd, posStart, posEnd);
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

        return compareChromosomes(Chromosome, other.Chromosome);
    }

    public static int compareChromosomes(final String chr1, final String chr2)
    {
        // we use chromosome rank such that chr1 is sorted before chr2
        int rank1 = HumanChromosome.chromosomeRank(chr1);
        int rank2 = HumanChromosome.chromosomeRank(chr2);

        if(rank1 == rank2)
        {
            // will occur for non-standard contigs, in which case revert to standard string comparison
            return chr1.compareTo(chr2);
        }

        return rank1 < rank2 ? -1 : 1;
    }

    public static List<ChrBaseRegion> loadChrBaseRegionList(final String filename)
    {
        Map<String,List<BaseRegion>> chrRegionMap = loadChrBaseRegions(filename);

        List<ChrBaseRegion> regions = Lists.newArrayList();

        for(Map.Entry<String,List<BaseRegion>> entry : chrRegionMap.entrySet())
        {
            entry.getValue().forEach(x -> regions.add(new ChrBaseRegion(entry.getKey(), x.start(), x.end())));
        }

        checkMergeOverlaps(regions, true);

        return regions;
    }

    public static Map<String,List<BaseRegion>> loadChrBaseRegions(final String filename)
    {
        if(filename == null)
            return Collections.emptyMap();

        boolean isBedFile = filename.endsWith(".bed") || filename.endsWith(".bed.gz");
        return loadChrBaseRegions(filename, isBedFile);
    }

    public static Map<String,List<BaseRegion>> loadChrBaseRegions(final String filename, boolean isBedFile)
    {
        Map<String,List<BaseRegion>> chrRegionsMap = Maps.newHashMap();

        if(filename == null)
            return chrRegionsMap;

        // accepts zipped / non-zipped, with and without headers, Chromosome/chromosome,PosStart/PositionStart etc
        try(BufferedReader fileReader = createBufferedReader(filename))
        {
            String delim = FileDelimiters.inferFileDelimiter(filename);

            int chrIndex = 0;
            int posStartIndex = 1;
            int posEndIndex = 2;

            String line = fileReader.readLine();

            boolean hasHeader = line.contains(FLD_CHROMOSOME);

            if(hasHeader)
            {
                Map<String,Integer> fieldIndexMap = FileReaderUtils.createFieldsIndexMap(line, delim);

                chrIndex = getChromosomeFieldIndex(fieldIndexMap);
                posStartIndex = getPositionStartFieldIndex(fieldIndexMap);
                posEndIndex = getPositionEndFieldIndex(fieldIndexMap);

                line = fileReader.readLine();
            }

            while(line != null)
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

                line = fileReader.readLine();
            }

            chrRegionsMap.values().forEach(x -> BaseRegion.checkMergeOverlaps(x));

            return chrRegionsMap;
        }
        catch(IOException e)
        {
            return null;
        }
    }

    public static void checkMergeOverlaps(final List<ChrBaseRegion> regions, boolean checkSorted)
    {
        if(checkSorted)
            Collections.sort(regions);

        // merge any adjacent regions
        int index = 0;
        while(index < regions.size() - 1)
        {
            ChrBaseRegion region = regions.get(index);
            ChrBaseRegion nextRegion = regions.get(index + 1);

            if(region.Chromosome.equals(nextRegion.Chromosome) && region.end() >= nextRegion.start() - 2)
            {
                region.setEnd(nextRegion.end());
                regions.remove(index + 1);
            }
            else
            {
                ++index;
            }
        }
    }

    public static final int INVALID_FIELD = -1;

    public static int getChromosomeFieldIndex(final Map<String,Integer> fieldIndexMap)
    {
        if(fieldIndexMap.containsKey(FLD_CHROMOSOME))
            return fieldIndexMap.get(FLD_CHROMOSOME);
        else if(fieldIndexMap.containsKey(FLD_CHROMOSOME.toLowerCase()))
            return fieldIndexMap.get(FLD_CHROMOSOME.toLowerCase());

        return INVALID_FIELD;
    }

    public static int getPositionFieldIndex(final Map<String,Integer> fieldIndexMap)
    {
        if(fieldIndexMap.containsKey(FLD_POSITION))
            return fieldIndexMap.get(FLD_POSITION);
        else if(fieldIndexMap.containsKey(FLD_POSITION.toLowerCase()))
            return fieldIndexMap.get(FLD_POSITION.toLowerCase());

        return INVALID_FIELD;
    }

    public static int getPositionStartFieldIndex(final Map<String,Integer> fieldIndexMap)
    {
        if(fieldIndexMap.containsKey(FLD_POSITION_START))
            return fieldIndexMap.get(FLD_POSITION_START);
        else if(fieldIndexMap.containsKey(FLD_POSITION_START.toLowerCase()))
            return fieldIndexMap.get(FLD_POSITION_START.toLowerCase());

        if(fieldIndexMap.containsKey(FLD_POS_START))
            return fieldIndexMap.get(FLD_POS_START);
        else if(fieldIndexMap.containsKey(FLD_POS_START.toLowerCase()))
            return fieldIndexMap.get(FLD_POS_START.toLowerCase());

        return INVALID_FIELD;
    }

    public static int getPositionEndFieldIndex(final Map<String,Integer> fieldIndexMap)
    {
        if(fieldIndexMap.containsKey(FLD_POSITION_END))
            return fieldIndexMap.get(FLD_POSITION_END);
        else if(fieldIndexMap.containsKey(FLD_POSITION_END.toLowerCase()))
            return fieldIndexMap.get(FLD_POSITION_END.toLowerCase());

        if(fieldIndexMap.containsKey(FLD_POS_END))
            return fieldIndexMap.get(FLD_POS_END);
        else if(fieldIndexMap.containsKey(FLD_POS_END.toLowerCase()))
            return fieldIndexMap.get(FLD_POS_END.toLowerCase());

        return INVALID_FIELD;
    }
}

