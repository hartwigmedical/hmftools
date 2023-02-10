package com.hartwig.hmftools.common.genome.unmappableregions;

import com.hartwig.hmftools.common.genome.chromosome.HumanChromosome;
import com.hartwig.hmftools.common.genome.refgenome.RefGenomeSource;
import htsjdk.samtools.reference.IndexedFastaSequenceFile;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.PrintWriter;
import java.util.ArrayList;

public class UnmappableRegionsCalculator
{
    private final RefGenomeSource mRefGenomeSource;
    private final String mOutputBedPath;
    private int mChunkSize = 1000000;

    // Constructors ================================
    public UnmappableRegionsCalculator(String refGenomeFastaPath, String outputBedPath) throws FileNotFoundException
    {
        mRefGenomeSource = new RefGenomeSource(new IndexedFastaSequenceFile(new File(refGenomeFastaPath)));
        mOutputBedPath = outputBedPath;
    }

    // Setters ================================
    public UnmappableRegionsCalculator setChunkSize(int chunkSize)
    {
        mChunkSize = chunkSize;
        return this;
    }

    // Methods ================================
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

    public Integer[][] chromosomeChunksStartEndPositions(String chromosome, int chunkSize) throws Exception
    {
        if(chunkSize < 1)
            throw new Exception("`chunkSize` must be >=1");

        int chromosomeLength = mRefGenomeSource.getChromosomeLength(chromosome);

        ArrayList<Integer> startPositions = new ArrayList<>();
        ArrayList<Integer> endPositions = new ArrayList<>();

        for(int i = 0; i < chromosomeLength; i += chunkSize)
        {
            startPositions.add(i); // zero based positions
            endPositions.add(i + chunkSize);
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

    public ArrayList<Integer[]> unmappablePositionsInChromosome(String chromosome, int chunkSize) throws Exception
    {
        Integer[][] chunksStartEndPositions = chromosomeChunksStartEndPositions(chromosome, chunkSize);

        ArrayList<Integer[]> unmappedStartEndPositions = new ArrayList<>();

        for(int i = 0; i < chunksStartEndPositions.length; i++)
        {
            Integer chunkStartPosition = chunksStartEndPositions[i][0];
            Integer chunkEndPosition = chunksStartEndPositions[i][1];

            String chunkDnaSequence = mRefGenomeSource.getBaseString(
                chromosome,
                chunkStartPosition + 1, // getBaseString() uses 1 based positions: start and end are inclusive
                chunkEndPosition
            );

            Integer[][] currentUnmappedPositions = nucleotideStretchStartEndPositions(chunkDnaSequence, 'N');

            for(Integer[] j : currentUnmappedPositions)
            {
                // Convert chunk position to genome position
                j[0] += chunkStartPosition;
                j[1] += chunkStartPosition;
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
            String chromosome = humanChromosome.toString();
            System.out.println("Calculating unmapped regions for chromosome: " + chromosome);

            ArrayList<Integer[]> unmappedStartEndPositions = unmappablePositionsInChromosome(chromosome, mChunkSize);
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
        String refGenomeFastaPath = "/home/lnguyen/Documents/BamMetricsTest/resources/HMFtools-Resources_ref_genome_37_Homo_sapiens.GRCh37.GATK.illumina.fasta";
        String outputBedPath = "/home/lnguyen/Documents/BamMetricsTest/output/GRCh37_unmapped_regions.bed";
        UnmappableRegionsCalculator unmappableRegionsCalculator = new UnmappableRegionsCalculator(refGenomeFastaPath, outputBedPath);

//        unmappableRegionsCalculator.nucleotideStretchStartEndPositions("NNNAAANNNNNNTTT",'N',false);
//        unmappableRegionsCalculator.chromosomeChunksStartEndPositions("21", 1000000, true);

//        ArrayList<Integer[]> startEndPositions = unmappableRegionsCalculator.unmappablePositionsInChromosome("21", 2000000);
//        unmappableRegionsCalculator.flattenStartEndPositions(startEndPositions);

        unmappableRegionsCalculator.run();
    }
}
