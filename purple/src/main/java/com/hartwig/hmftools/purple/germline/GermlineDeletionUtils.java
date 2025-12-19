package com.hartwig.hmftools.purple.germline;

import static com.hartwig.hmftools.common.purple.GermlineStatus.HET_DELETION;
import static com.hartwig.hmftools.common.purple.GermlineStatus.HOM_DELETION;
import static com.hartwig.hmftools.common.region.BaseRegion.positionsOverlap;
import static com.hartwig.hmftools.purple.PurpleUtils.PPL_LOGGER;

import java.io.IOException;
import java.util.List;
import java.util.Map;
import java.util.stream.Collectors;

import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import com.hartwig.hmftools.common.gene.ExonData;
import com.hartwig.hmftools.common.gene.GeneData;
import com.hartwig.hmftools.common.gene.TranscriptData;
import com.hartwig.hmftools.common.purple.PurpleCopyNumber;
import com.hartwig.hmftools.common.purple.PurpleCopyNumberFile;
import com.hartwig.hmftools.common.purple.PurpleSegment;
import com.hartwig.hmftools.purple.region.ObservedRegion;

/**
 * Utility class for common operations related to germline deletions.
 * Shared functionality between SampleGermlineGeneTask and GermlineDeletions.
 */
public final class GermlineDeletionUtils
{
    private GermlineDeletionUtils()
    {
        // Prevent instantiation of utility class
    }

    /**
     * Loads Purple data (copy numbers and fitted regions) for a sample.
     *
     * @param sampleId The sample ID
     * @param purpleDir The Purple data directory
     * @param copyNumberMap Map to populate with copy number data
     * @param fittedRegionMap Map to populate with fitted region data
     */
    public static void loadPurpleDataForGermline(
            final String sampleId, final String purpleDir,
            final Map<String, List<PurpleCopyNumber>> copyNumberMap,
            final Map<String, List<ObservedRegion>> fittedRegionMap)
    {
        String samplePurpleDir = purpleDir.contains("*") ? purpleDir.replaceAll("\\*", sampleId) : purpleDir;

        try
        {
            List<PurpleCopyNumber> allCopyNumbers = PurpleCopyNumberFile.read(
                    PurpleCopyNumberFile.generateFilenameForReading(samplePurpleDir, sampleId));

            List<PurpleSegment> segments = PurpleSegment.read(PurpleSegment.generateFilename(samplePurpleDir, sampleId)).stream()
                    .filter(x -> x.GermlineState == HET_DELETION || x.GermlineState == HOM_DELETION)
                    .collect(Collectors.toList());

            for(PurpleSegment segment : segments)
            {
                List<ObservedRegion> regions = fittedRegionMap.get(segment.Chromosome);

                if(regions == null)
                {
                    regions = Lists.newArrayList();
                    fittedRegionMap.put(segment.Chromosome, regions);

                    copyNumberMap.put(
                            segment.Chromosome,
                            allCopyNumbers.stream().filter(x -> x.chromosome().equals(segment.Chromosome)).collect(Collectors.toList()));
                }

                regions.add(ObservedRegion.fromSegment(segment));
            }

            PPL_LOGGER.debug("sample({}) read {} het-hom deletion regions",
                    sampleId, fittedRegionMap.values().stream().mapToInt(x -> x.size()).sum());
        }
        catch(IOException e)
        {
            PPL_LOGGER.error("sample({}) failed to load purple files form {}: {}", sampleId, samplePurpleDir, e.toString());
        }
    }

    /**
     * Finds a matching copy number for a region.
     *
     * @param region The region to find a matching copy number for
     * @param copyNumbers List of copy numbers to search
     * @return The matching copy number, or null if none found
     */
    public static PurpleCopyNumber findMatchingCopyNumber(
            final ObservedRegion region, final List<PurpleCopyNumber> copyNumbers)
    {
        if(copyNumbers == null || copyNumbers.isEmpty())
        {
            return null;
        }

        String chromosome = region.chromosome();
        int regionStart = region.start();
        int regionEnd = region.end();

        // Find all overlapping copy numbers
        List<PurpleCopyNumber> overlappingCopyNumbers = copyNumbers.stream()
                .filter(x -> x.chromosome().equals(chromosome) && positionsOverlap(x.start(), x.end(), regionStart, regionEnd))
                .collect(Collectors.toList());

        if(overlappingCopyNumbers.isEmpty())
        {
            return null;
        }

        // If there's only one overlapping copy number, return it
        if(overlappingCopyNumbers.size() == 1)
        {
            return overlappingCopyNumbers.get(0);
        }

        // If there are multiple overlapping copy numbers, prefer the one with the lowest major allele copy number
        // This matches the behavior in the original GermlineDeletions class
        return overlappingCopyNumbers.stream()
                .min((a, b) -> Double.compare(a.majorAlleleCopyNumber(), b.majorAlleleCopyNumber()))
                .orElse(overlappingCopyNumbers.get(0));
    }

    /**
     * Finds genes that overlap with a region.
     *
     * @param chromosome The chromosome
     * @param regionStart The region start position
     * @param regionEnd The region end position
     * @param buffer Buffer to add to region boundaries
     * @param geneDataList List of genes to search
     * @return List of overlapping genes
     */
    public static List<GeneData> findOverlappingGenes(
            final String chromosome, final int regionStart, final int regionEnd, final int buffer,
            final List<GeneData> geneDataList)
    {
        if(geneDataList == null || geneDataList.isEmpty())
        {
            return Lists.newArrayList();
        }

        int regionLowerPos = regionStart - buffer;
        int regionHighPos = regionEnd + buffer;

        return geneDataList.stream()
                .filter(x -> positionsOverlap(x.GeneStart, x.GeneEnd, regionLowerPos, regionHighPos))
                .collect(Collectors.toList());
    }

    /**
     * Finds exons that overlap with a region.
     *
     * @param transData The transcript data containing exons
     * @param regionStart The region start position
     * @param regionEnd The region end position
     * @param buffer Buffer to add to region boundaries
     * @return List of overlapping exons
     */
    public static List<ExonData> findOverlappingExons(
            final TranscriptData transData, final int regionStart, final int regionEnd, final int buffer)
    {
        if(transData == null || transData.exons().isEmpty())
        {
            return Lists.newArrayList();
        }

        int regionLowerPos = regionStart - buffer;
        int regionHighPos = regionEnd + buffer;

        return transData.exons().stream()
                .filter(x -> positionsOverlap(x.Start, x.End, regionLowerPos, regionHighPos))
                .collect(Collectors.toList());
    }
}
