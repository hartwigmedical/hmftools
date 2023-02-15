package com.hartwig.hmftools.bamtools.unmappableregions;

import com.hartwig.hmftools.common.genome.chromosome.HumanChromosome;
import com.hartwig.hmftools.common.genome.refgenome.RefGenomeCoordinates;
import com.hartwig.hmftools.common.genome.refgenome.RefGenomeSource;

import static com.hartwig.hmftools.bamtools.common.CommonUtils.BT_LOGGER;
import static com.hartwig.hmftools.bamtools.common.CommonUtils.PARTITION_SIZE;
import static com.hartwig.hmftools.common.genome.refgenome.RefGenomeFunctions.stripChrPrefix;
import static com.hartwig.hmftools.common.genome.refgenome.RefGenomeSource.REF_GENOME;
import static com.hartwig.hmftools.common.genome.refgenome.RefGenomeSource.addRefGenomeConfig;

import com.hartwig.hmftools.common.genome.refgenome.RefGenomeVersion;
import htsjdk.samtools.reference.IndexedFastaSequenceFile;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.List;

import static com.hartwig.hmftools.bamtools.common.CommonUtils.DEFAULT_CHR_PARTITION_SIZE;
import static com.hartwig.hmftools.common.genome.refgenome.RefGenomeVersion.*;
import static com.hartwig.hmftools.common.utils.FileWriterUtils.*;

import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.CommandLineParser;
import org.apache.commons.cli.DefaultParser;
import org.apache.commons.cli.Options;

public class UnmappableRegionsCalculator
{
    private final RefGenomeSource mRefGenomeSource;
    private final IndexedFastaSequenceFile mIndexedFastaSequenceFile;

    private final String mOutputBedPath;
    private int mPartitionSize = DEFAULT_CHR_PARTITION_SIZE; // 1000000
    private RefGenomeVersion mRefGenomeVersion = V37;

    public UnmappableRegionsCalculator(String refGenomeFastaPath, String outputBedPath) throws FileNotFoundException
    {
        mIndexedFastaSequenceFile = new IndexedFastaSequenceFile(new File(refGenomeFastaPath));
        mRefGenomeSource = new RefGenomeSource(mIndexedFastaSequenceFile);
        mOutputBedPath = outputBedPath;
    }

    public void setPartitionSize(int partitionSize)
    {
        mPartitionSize = partitionSize;
    }

    public void setRefGenomeVersion(RefGenomeVersion refGenomeVersion)
    {
        mRefGenomeVersion = refGenomeVersion;
    }

    public Integer[][] nucleotideStretchStartEndPositions(String nucleotideSequence, char targetNucleotide) throws Exception
    {
//        String nucleotideSequence = mRefGenomeSource.getBaseString("21",1,10);
//        String nucleotideSequence = "NNNNN";
//        System.out.println(nucleotideSequence);
//        System.out.println(nucleotideSequence.length());
//        char targetNucleotide = 'N';

        if(nucleotideSequence.length() > Integer.MAX_VALUE)
            throw new Exception("`nucleotideSequence` exceeds max length of " + Integer.MAX_VALUE);

        ArrayList<Integer> startPositions = new ArrayList<>();
        ArrayList<Integer> endPositions = new ArrayList<>();

        boolean isNucleotideStretch = false;
        for(int i = 0; i < nucleotideSequence.length(); i++)
        {
            char nucleotide = nucleotideSequence.charAt(i); // charAt accepts int and not long

            if(nucleotide == targetNucleotide)
            {
                if(!isNucleotideStretch)
                {
                    isNucleotideStretch = true;
                    startPositions.add(i); // zero based positions
                }
            } else {
                if(isNucleotideStretch)
                {
                    isNucleotideStretch = false;
                    endPositions.add(i);
                }
            }
        }

        // The last end position should be the length of `nucleotideSequence` if the tail of the sequence is still a
        // nucleotide stretch
        if(isNucleotideStretch)
            endPositions.add(nucleotideSequence.length());

        // Output array of pairs of start/end positions
        Integer[][] startEndPositions = new Integer[startPositions.size()][2];
        for(int i = 0; i < startEndPositions.length; i++)
        {
            startEndPositions[i][0] = startPositions.get(i);
            startEndPositions[i][1] = endPositions.get(i);
        }

//        // For debug
//        System.out.println(nucleotideSequence);
//        for(Integer[] i : startEndPositions)
//            System.out.println(Arrays.toString(i));

        return startEndPositions;
    }

    public Integer[][] chromosomePartitionsStartEndPositions(String chromosome, int partitionSize) throws Exception
    {
        if(partitionSize < 1)
            throw new Exception("`partitionSize` must be >=1");

        RefGenomeCoordinates refGenomeCoords = mRefGenomeVersion == V37 ? RefGenomeCoordinates.COORDS_37 : RefGenomeCoordinates.COORDS_38;
        int chromosomeLength = refGenomeCoords.length(stripChrPrefix(chromosome));
//        int chromosomeLength = mRefGenomeSource.getChromosomeLength(chromosome);

        ArrayList<Integer> startPositions = new ArrayList<>();
        ArrayList<Integer> endPositions = new ArrayList<>();

        for(int i = 0; i < chromosomeLength; i += partitionSize)
        {
            startPositions.add(i); // zero based positions
            endPositions.add(i + partitionSize);
        }

        // Set last end position to be the chrom length
        endPositions.set(endPositions.size() - 1, chromosomeLength);

        // Output array of pairs of start/end positions
        Integer[][] startEndPositions = new Integer[startPositions.size()][2];
        for(int i = 0; i < startEndPositions.length; i++)
        {
            startEndPositions[i][0] = startPositions.get(i);
            startEndPositions[i][1] = endPositions.get(i);
        }

//        // For debug
//        for(Integer[] i : startEndPositions)
//            System.out.println(Arrays.toString(i));

        return startEndPositions;
    }

    public ArrayList<Integer[]> unmappablePositionsInChromosome(HumanChromosome humanChromosome, int partitionSize) throws Exception
    {
        String chromosome = mRefGenomeVersion.versionedChromosome(humanChromosome.toString());
        ArrayList<Integer[]> output = unmappablePositionsInChromosome(chromosome, partitionSize);
        return output;
    }

    public ArrayList<Integer[]> unmappablePositionsInChromosome(String chromosome, int partitionSize) throws Exception
    {
        Integer[][] partitionsStartEndPositions = chromosomePartitionsStartEndPositions(chromosome, partitionSize);
        ArrayList<Integer[]> unmappedStartEndPositions = new ArrayList<>();

        for(int i = 0; i < partitionsStartEndPositions.length; i++)
        {
            Integer partitionStartPosition = partitionsStartEndPositions[i][0];
            Integer partitionEndPosition = partitionsStartEndPositions[i][1];

            String partitionDnaSequence = mRefGenomeSource.getBaseString(
                chromosome,
                partitionStartPosition + 1, // getBaseString() uses 1 based positions: start and end are inclusive
                partitionEndPosition
            );

            Integer[][] currentUnmappedPositions = nucleotideStretchStartEndPositions(partitionDnaSequence, 'N');

            for(Integer[] j : currentUnmappedPositions)
            {
                // Convert partition position to genome position
                j[0] += partitionStartPosition;
                j[1] += partitionStartPosition;
                //System.out.println(Arrays.toString(j));

                unmappedStartEndPositions.add(j);
            }
        }

//        // Debug
//        for(Integer[] i : unmappedStartEndPositions)
//        {
//            System.out.println(Arrays.toString(i));
//        }
//        System.out.println();

//        int start = 9595548;
//        int end = 9645548;
////        System.out.println(mRefGenomeSource.getBaseString(chromosome, start+1, end)); // Must be all Ns
//        System.out.println(mRefGenomeSource.getBaseString(chromosome, start+1-1, end+1)); // All Ns with one nucleotide on either side

        return unmappedStartEndPositions;
    }

    public ArrayList<Integer[]> flattenStartEndPositions(ArrayList<Integer[]> startEndPositions)
    {
        // Group indexes of `startEndPositions` that are contiguous
        ArrayList<ArrayList<Integer>> contiguousRegionIndexGroups = new ArrayList<>();
        ArrayList<Integer> currentContiguousRegionIndexes = new ArrayList<>();
        for(int i = 0; i < startEndPositions.size(); i++)
        {
            if(i==0)
            {
                currentContiguousRegionIndexes.add(i);
                continue;
            }

            Integer previousEnd = startEndPositions.get(i-1)[1];
            Integer currentStart = startEndPositions.get(i)[0];

            if(!previousEnd.equals(currentStart))
            {
                contiguousRegionIndexGroups.add(currentContiguousRegionIndexes);
                currentContiguousRegionIndexes = new ArrayList<>();
            }

            // Implied: if(previousEnd.equals(currentStart))
            currentContiguousRegionIndexes.add(i);
        }

//        // Debug
//        for(ArrayList<Integer> i : contiguousRegionIndexGroups)
//            System.out.println(i.toString());
//        System.out.println();

        // Get start and end positions of contiguousRegions
        ArrayList<Integer[]> startEndPositionsFlattened = new ArrayList<>();
        for(ArrayList<Integer> indexGroup : contiguousRegionIndexGroups)
        {
            Integer groupFirstIndex = indexGroup.get(0);
            Integer groupLastIndex;
            Integer[] groupStartEndPositions = new Integer[2];

            if(indexGroup.size() == 1)
            {
                groupStartEndPositions = startEndPositions.get( groupFirstIndex ); // get both start and end
            } else {
                groupLastIndex = indexGroup.get( indexGroup.size()-1 );
                groupStartEndPositions[0] = startEndPositions.get(groupFirstIndex)[0]; // start
                groupStartEndPositions[1] = startEndPositions.get(groupLastIndex)[1]; // end
            }

            startEndPositionsFlattened.add(groupStartEndPositions);
        }

//        // Debug
//        for(Integer[] i : startEndPositionsFlattened)
//            System.out.println(Arrays.toString(i));
//        System.out.println();

        return startEndPositionsFlattened;
    }

    public void run() throws Exception
    {
        PrintWriter writer = new PrintWriter(mOutputBedPath);
        writer.println("chrom\tstart\tend");

        for(HumanChromosome humanChromosome : HumanChromosome.values())
        {
//            String chromosome = humanChromosome.toString();
            String chromosome = mRefGenomeVersion.versionedChromosome(humanChromosome.toString());

            BT_LOGGER.info("Calculating N-nucleotide regions for chromosome: " + chromosome);

            ArrayList<Integer[]> unmappedStartEndPositions = unmappablePositionsInChromosome(humanChromosome, mPartitionSize);
            ArrayList<Integer[]> unmappedStartEndPositionsFlattened = flattenStartEndPositions(unmappedStartEndPositions);

            for(Integer[] startEndPosition : unmappedStartEndPositionsFlattened)
            {
                writer.println(chromosome + "\t" + startEndPosition[0] + "\t" + startEndPosition[1]);
            }
        }

        writer.close();
    }

    // Entry point ================================
    public static void main(String[] args) throws Exception
    {
        // Initialize options
        Options options = new Options();
        CommandLineParser parser = new DefaultParser();

        // Add options
        addRefGenomeConfig(options); // -ref_genome, -ref_genome_version
        addOutputDir(options); // -output_dir,
        options.addOption(PARTITION_SIZE, true, "Partition size, default: " + DEFAULT_CHR_PARTITION_SIZE);
        CommandLine cmd = parser.parse(options, args);

        // Apply commandline arguments
        String outputBedPath = cmd.getOptionValue(OUTPUT_DIR) + "genome_unmappable_regions." + cmd.getOptionValue(REF_GENOME_VERSION) + ".bed";

        UnmappableRegionsCalculator unmappableRegionsCalculator = new UnmappableRegionsCalculator(
            cmd.getOptionValue(REF_GENOME),
            outputBedPath
        );

        unmappableRegionsCalculator.setRefGenomeVersion(
            RefGenomeVersion.from(cmd.getOptionValue(REF_GENOME_VERSION))
        );

        int partitionSize = Integer.parseInt(cmd.getOptionValue(PARTITION_SIZE, String.valueOf(DEFAULT_CHR_PARTITION_SIZE)));
        unmappableRegionsCalculator.setPartitionSize(partitionSize);

        //
        unmappableRegionsCalculator.run();
    }
}
