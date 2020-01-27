package com.hartwig.hmftools.svtools.rna_expression;

import static java.lang.Math.max;
import static java.lang.Math.min;

import java.io.File;
import java.io.IOException;
import java.util.List;
import java.util.stream.Collectors;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.genome.region.GenomeRegion;
import com.hartwig.hmftools.common.genome.region.GenomeRegions;
import com.hartwig.hmftools.common.variant.hotspot.SAMSlicer;

import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.Options;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.jetbrains.annotations.NotNull;

import htsjdk.samtools.CigarOperator;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;
import htsjdk.samtools.reference.IndexedFastaSequenceFile;

public class RnaBamReader
{
    private IndexedFastaSequenceFile mIndexedFastaSequenceFile;
    private final File mRefGenomeFile;
    private final String mBamFile;
    private final SamReader mSamReader;

    private final List<SAMRecord> mBamRecords;

    private static final int DEFAULT_MIN_BASE_QUALITY = 13;
    private static final int DEFAULT_MIN_MAPPING_QUALITY = 1;

    private static final Logger LOGGER = LogManager.getLogger(RnaBamReader.class);

    // config
    public static final String REF_GENOME = "ref_genome";
    public static final String BAM_FILE = "bam_file";


    public RnaBamReader(final CommandLine cmd)
    {
        mBamFile = cmd.getOptionValue(BAM_FILE);

        final String refGenomeFile = cmd.getOptionValue(REF_GENOME);
        mRefGenomeFile = new File(refGenomeFile);

        mBamRecords = Lists.newArrayList();

        mSamReader = SamReaderFactory.makeDefault().referenceSequence(mRefGenomeFile).open(new File(mBamFile));

        try
        {
            LOGGER.debug("Loading indexed fasta reference file");
            mIndexedFastaSequenceFile = new IndexedFastaSequenceFile(new File(refGenomeFile));
        }
        catch (IOException e)
        {
            LOGGER.error("Reference file loading failed");
            return;
        }

    }

    public static boolean validConfig(final CommandLine cmd)
    {
        return cmd.hasOption(REF_GENOME) && cmd.hasOption(BAM_FILE);
    }

    public static void addCommandLineOptions(final Options options)
    {
        options.addOption(REF_GENOME, true, "Ref genome file location");
        options.getOption(REF_GENOME).setRequired(true);
        options.addOption(BAM_FILE, true, "RNA BAM file location");
        options.getOption(BAM_FILE).setRequired(true);
    }

    public void readBamCounts(final GenomeRegion genomeRegion)
    {
        mBamRecords.clear();

        SAMSlicer samSlicer = new SAMSlicer(DEFAULT_MIN_MAPPING_QUALITY, Lists.newArrayList(genomeRegion));
        samSlicer.slice(mSamReader, this::processSamRecord);
    }

    private void processSamRecord(@NotNull final SAMRecord record)
    {
        mBamRecords.add(record);

        // selector.select(asRegion(record), bafEvidence -> bafFactory.addEvidence(bafEvidence, record));
    }

    public void analyseReads(final GeneReadData geneReadData)
    {
        // regionReadData.getBamRecords().addAll(mBamRecords);

        // cache reference bases for comparison with read bases
        for(RegionReadData region : geneReadData.getRegionReadData())
        {
            final String regionRefBases = mIndexedFastaSequenceFile.getSubsequenceAt(
                    region.chromosome(), region.start(), region.end()).getBaseString();

            region.setRefBases(regionRefBases);
        }

        // for each record find all exons with an overlap
        // skip records if either end isn't in one of the exons for this gene

        for(final SAMRecord record : mBamRecords)
        {
            // skip records if either end isn't in one of the genic regions
            boolean startMatched = geneReadData.getRegionReadData().stream()
                    .anyMatch(x -> record.getStart() >= x.Region.start() && record.getStart() <= x.Region.end());

            boolean endMatched = geneReadData.getRegionReadData().stream()
                    .anyMatch(x -> record.getEnd() >= x.Region.start() && record.getEnd() <= x.Region.end());

            if(!startMatched || !endMatched)
                continue;

            // the read covers a part of the exon
            List<RegionReadData> containingRegions = geneReadData.getRegionReadData().stream()
                    .filter(x -> x.Region.start() <= record.getStart() && x.Region.end() >= record.getEnd())
                    .collect(Collectors.toList());

            // the read covers all of the exon
            List<RegionReadData> containedRegions = geneReadData.getRegionReadData().stream()
                    .filter(x -> !containingRegions.contains(x))
                    .filter(x -> x.Region.start() >= record.getStart() && x.Region.end() <= record.getEnd())
                    .collect(Collectors.toList());

            // the read covers part of the exon and crosses its boundary
            List<RegionReadData> overlappingRegions = geneReadData.getRegionReadData().stream()
                    .filter(x -> !containedRegions.contains(x) && !containingRegions.contains(x))
                    .filter(x -> overlaps(x.Region, record))
                    .collect(Collectors.toList());

            if(containedRegions.isEmpty() && containingRegions.isEmpty() && overlappingRegions.isEmpty())
                continue;

            // record the links between these exons
            List<RegionReadData> allRegions = Lists.newArrayList();
            allRegions.addAll(containedRegions);
            allRegions.addAll(containingRegions);
            allRegions.addAll(overlappingRegions);

            if(allRegions.size() > 1)
            {
                for(int i = 0; i < allRegions.size() - 1; ++i)
                {
                    RegionReadData region1 = allRegions.get(i);

                    for(int j = i + 1; j < allRegions.size(); ++j)
                    {
                        RegionReadData region2 = allRegions.get(j);
                        region1.addLinkedRegion(region2);
                        region2.addLinkedRegion(region1);
                    }
                }
            }

            // now check for matching bases in the read vs the reference for the overlapping sections
            // final String recordBaseStr = record.getReadString();

            allRegions.forEach(x -> setMatchingBases(x, record));
        }

        // summary stats


        //LOGGER.debug("region({}:{}->{}): reads({}) overlaps({}) refBaseMatches({})",
        //        region.chromosome(), region.start(), region.end(), mBamRecords.size(), overlapMatchCount, refBaseMatchCount);
    }

    private boolean overlaps(final GenomeRegion region, final SAMRecord record)
    {
        return !(record.getStart() > region.end() || record.getEnd() < region.start());
    }

    private void setMatchingBases(final RegionReadData region, final SAMRecord record)
    {
        long overlapStart = max(region.start(), record.getStart());
        long overlapEnd = min(region.end(), record.getEnd());
        int overlapLength = (int)(overlapEnd - overlapStart + 1);

        if(overlapLength < 5)
            return;

        // compare the bases at this location
        int recordOverlapStart = (int)(overlapStart - record.getStart());

        // factor in soft-clipping
        if(record.getCigar().getFirstCigarElement().getOperator() == CigarOperator.S)
        {
            recordOverlapStart += record.getCigar().getFirstCigarElement().getLength();
        }

        int recordOverlapEnd = min(recordOverlapStart + overlapLength, record.getReadString().length() - 1);

        if(recordOverlapStart < 0 || recordOverlapStart >= record.getReadString().length())
        {
            recordOverlapStart = 0;
        }

        final String recordBases = record.getReadString().substring(recordOverlapStart, recordOverlapEnd);

        int regionOverlapStart = (int)(overlapStart - region.start());
        int regionOverlapEnd = regionOverlapStart + overlapLength;
        final String regionBases = region.refBases().substring(regionOverlapStart, regionOverlapEnd);

        if(!recordBases.equals(regionBases))
        {
            // prform a manual comparison
            int matchedBases = findStringOverlaps(regionBases, recordBases);

            if(matchedBases < MIN_BASE_MATCH_PERC * overlapLength)
            {
                LOGGER.trace("region({}) has base-mismatch, overlap({}) matched({}) cigar({})",
                        region, overlapLength, matchedBases, record.getCigar());

                LOGGER.trace("regionBases: pos({} -> {}) {}",
                        regionOverlapStart, regionOverlapEnd, regionBases);
                LOGGER.trace("recordBases: pos({} -> {}) {}",
                        recordOverlapStart, recordOverlapEnd, recordBases);
                return;
            }
        }

        int[] matchedBases = region.refBasesMatched();

        for(long i = overlapStart; i <= overlapEnd; ++i)
        {
            int baseIndex = (int)(i - region.start());
            ++matchedBases[baseIndex];
        }

        region.addMatchedRead();
    }

    private static final double MIN_BASE_MATCH_PERC = 0.9;

    private int findStringOverlaps(final String str1, final String str2)
    {
        if(str1.length() == 0 || str2.length() == 0)
            return 0;

        int matched = 0;
        int i = 0;
        int j = 0;
        int mismatchIndex = -1;

        // first compare bases at same indices, making note of the first difference if there is one
        while(i < str1.length() && j < str2.length())
        {
            if (str1.charAt(i) == str2.charAt(j))
                ++matched;
            else if(mismatchIndex == -1)
                mismatchIndex = i;

            ++i;
            ++j;
        }

        if(matched > MIN_BASE_MATCH_PERC * min(str1.length(), str2.length()))
            return matched;

        i = j = mismatchIndex;
        matched = mismatchIndex;

        while(i < str1.length() && j < str2.length())
        {
            if(str1.charAt(i) == str2.charAt(j))
            {
                ++i;
                ++j;
                ++matched;
                continue;
            }

            // search ahead in each string in turn for the next short matching sequence
            int startI = i;
            boolean seqFound = false;
            for(; i < str1.length() - 2 && j < str2.length() - 2; ++i)
            {
                if(str1.charAt(i) == str2.charAt(j) && str1.charAt(i+1) == str2.charAt(j+1) && str1.charAt(i+2) == str2.charAt(j+2))
                {
                    seqFound = true;
                    break;
                }
            }

            if(seqFound)
                continue;

            i = startI;

            for(; i < str1.length() - 2 && j < str2.length() - 2; ++j)
            {
                if(str1.charAt(i) == str2.charAt(j) && str1.charAt(i+1) == str2.charAt(j+1) && str1.charAt(i+2) == str2.charAt(j+2))
                {
                    seqFound = true;
                    break;
                }
            }

            if(!seqFound)
                break;
        }

        return matched;
    }

    public void testReadBamCounts()
    {
        LOGGER.debug("reading BAM file: {}", mBamFile);

        // slicing
        List<GenomeRegion> genomeRegions = Lists.newArrayList();

        String chromosome = "21";
        long posStart = 42870046;
        long posEnd = 42870116;

        genomeRegions.add(GenomeRegions.create(chromosome, posStart, posEnd));

        SAMSlicer samSlicer = new SAMSlicer(DEFAULT_MIN_MAPPING_QUALITY, genomeRegions);

        samSlicer.slice(mSamReader, this::processSamRecord);

        // mTumorReader.


    }



            /*
        for(BachelorGermlineVariant variant : bachRecords)
        {
            VariantHotspot variantHotspot = ImmutableVariantHotspotImpl.builder()
                    .chromosome(variant.Chromosome)
                    .position(variant.Position)
                    .ref(variant.Ref)
                    .alt(variant.Alts)
                    .build();

            allHotspots.add(variantHotspot);
        }

        final VariantHotspotEvidenceFactory hotspotEvidenceFactory = new VariantHotspotEvidenceFactory(DEFAULT_MIN_MAPPING_QUALITY, DEFAULT_MIN_BASE_QUALITY, allHotspots);
        final List<VariantHotspotEvidence> tumorEvidence = hotspotEvidenceFactory.evidence(mIndexedFastaSequenceFile, mSamReader);

        if(tumorEvidence.size() != bachRecords.size())
        {
            LOGGER.error("Incomplete BAM evidence read: evidenceCount({}) vs bachRecords({})", tumorEvidence.size(), bachRecords.size());
            return;
        }


        for(BachelorGermlineVariant variant : bachRecords)
        {
            for(VariantHotspotEvidence evidence : tumorEvidence)
            {
                if(evidence.chromosome().equals(variant.Chromosome) && evidence.position() == variant.Position
                        && evidence.ref().equals(variant.Ref) && evidence.alt().equals(variant.Alts))
                {
                    variant.setTumorData(evidence.altSupport(), evidence.readDepth());

                    LOGGER.debug("chr({}) position({}) matched, counts(ref={} alt={} depth={})",
                            variant.Chromosome, variant.Position,
                            variant.getTumorRefCount(), variant.getTumorAltCount(), variant.getTumorReadDepth());

                    break;
                }
            }
        }

     */
}
