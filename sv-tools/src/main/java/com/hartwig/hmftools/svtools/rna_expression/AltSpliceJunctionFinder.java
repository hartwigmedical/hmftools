package com.hartwig.hmftools.svtools.rna_expression;

import static java.lang.Math.abs;
import static java.lang.Math.min;

import static com.hartwig.hmftools.common.utils.io.FileWriterUtils.createBufferedWriter;
import static com.hartwig.hmftools.linx.types.SvVarData.SE_END;
import static com.hartwig.hmftools.linx.types.SvVarData.SE_PAIR;
import static com.hartwig.hmftools.linx.types.SvVarData.SE_START;
import static com.hartwig.hmftools.svtools.rna_expression.AltSpliceJunction.CONTEXT_EXONIC;
import static com.hartwig.hmftools.svtools.rna_expression.AltSpliceJunction.CONTEXT_INTRONIC;
import static com.hartwig.hmftools.svtools.rna_expression.AltSpliceJunction.CONTEXT_MIXED;
import static com.hartwig.hmftools.svtools.rna_expression.AltSpliceJunction.CONTEXT_SJ;
import static com.hartwig.hmftools.svtools.rna_expression.AltSpliceJunctionType.CRYPTIC_SPLICE_SITE;
import static com.hartwig.hmftools.svtools.rna_expression.AltSpliceJunctionType.NEW_EXONS;
import static com.hartwig.hmftools.svtools.rna_expression.AltSpliceJunctionType.SKIPPED_EXONS;
import static com.hartwig.hmftools.svtools.rna_expression.GeneBamReader.DEFAULT_MIN_MAPPING_QUALITY;
import static com.hartwig.hmftools.svtools.rna_expression.RnaExpConfig.RE_LOGGER;
import static com.hartwig.hmftools.svtools.rna_expression.RnaExpUtils.positionWithin;
import static com.hartwig.hmftools.svtools.rna_expression.RnaExpUtils.positionsOverlap;
import static com.hartwig.hmftools.svtools.rna_expression.TransMatchType.SPLICE_JUNCTION;

import java.io.BufferedWriter;
import java.io.IOException;
import java.util.List;
import java.util.Map;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.genome.region.GenomeRegion;
import com.hartwig.hmftools.common.genome.region.GenomeRegions;
import com.hartwig.hmftools.common.variant.hotspot.SAMSlicer;
import com.hartwig.hmftools.common.variant.structural.annotation.TranscriptData;

import org.jetbrains.annotations.NotNull;

import htsjdk.samtools.CigarOperator;
import htsjdk.samtools.QueryInterval;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SamReader;

public class AltSpliceJunctionFinder
{
    private final RnaExpConfig mConfig;
    private final SamReader mSamReader;

    private final List<AltSpliceJunction> mAltSpliceJunctions;
    private final BufferedWriter mWriter;

    private GeneReadData mGene;
    private long mCurrentSjPosition;

    public AltSpliceJunctionFinder(final RnaExpConfig config, final SamReader samReader, final BufferedWriter writer)
    {
        mConfig = config;
        mSamReader = samReader;
        mAltSpliceJunctions = Lists.newArrayList();
        mWriter = writer;
        mCurrentSjPosition = 0;
        mGene = null;
    }

    public List<AltSpliceJunction> getAltSpliceJunctions() { return mAltSpliceJunctions; }

    public void setGeneData(final GeneReadData gene)
    {
        mGene = gene;
        mAltSpliceJunctions.clear();
    }

    public static boolean isCandidate(final ReadRecord read)
    {
        if(!read.Cigar.containsOperator(CigarOperator.N))
            return false;

        if(read.getTranscriptClassifications().values().contains(SPLICE_JUNCTION))
            return false;

        if(read.getMappedRegionCoords().size() == 1)
            return false;

        return true;
    }

    public boolean junctionMatchesGene(final long[] spliceJunction, final List<Integer> transIds)
    {
        for(TranscriptData transData : mGene.getTranscripts())
        {
            if(!transIds.contains(transData.TransId))
                continue;

            for(int i = 0; i < transData.exons().size() - 1; ++i)
            {
                if(transData.exons().get(i).ExonEnd == spliceJunction[SE_START] && transData.exons().get(i+1).ExonStart == spliceJunction[SE_END])
                    return true;
            }
        }

        // also check any overlapping gene's exonic regions
        if(mGene.getOtherGeneExonicRegions().stream()
                .anyMatch(x -> positionsOverlap(spliceJunction[SE_START], spliceJunction[SE_END], x[SE_START], x[SE_END])))
        {
            return true;
        }

        return false;
    }

    private AltSpliceJunction createFromRead(final ReadRecord read, final List<Integer> relatedTransIds)
    {
        // related transcripts will any of those where either read covers one or more of its exons
        long[] spliceJunction = new long[SE_PAIR];

        // find the novel splice junction, and all associated transcripts
        final List<long[]> mappedCoords = read.getMappedRegionCoords();
        spliceJunction[SE_START] = mappedCoords.get(0)[SE_END];
        spliceJunction[SE_END] = mappedCoords.get(1)[SE_START];

        List<RegionReadData> sjStartRegions = Lists.newArrayList(); // transcript regions with an exon matching the start of the alt SJ
        List<RegionReadData> sjEndRegions = Lists.newArrayList();

        String startContext = "";
        String endContext = "";

        for (Map.Entry<RegionReadData, RegionMatchType> entry : read.getMappedRegions().entrySet())
        {
            final RegionReadData region = entry.getKey();
            RegionMatchType matchType = entry.getValue();

            String regionStartContext = CONTEXT_INTRONIC;
            String regionEndContext = CONTEXT_INTRONIC;

            if (matchType == RegionMatchType.EXON_BOUNDARY || matchType == RegionMatchType.EXON_MATCH)
            {
                if (region.end() == spliceJunction[SE_START])
                {
                    startContext = CONTEXT_SJ;
                    sjStartRegions.add(region);
                }

                if (region.start() == spliceJunction[SE_END])
                {
                    endContext = CONTEXT_SJ;
                    sjEndRegions.add(region);
                }
            }
            else if(matchType == RegionMatchType.WITHIN_EXON)
            {
                if(positionWithin(spliceJunction[SE_START], region.start(), region.end()))
                {
                    regionStartContext = CONTEXT_EXONIC;
                }

                if(positionWithin(spliceJunction[SE_END], region.start(), region.end()))
                {
                    regionEndContext = CONTEXT_EXONIC;
                }
            }

            if(startContext != CONTEXT_SJ)
                startContext = startContext == "" ? regionStartContext : CONTEXT_MIXED;

            if(endContext != CONTEXT_SJ)
                endContext = endContext == "" ? regionEndContext : CONTEXT_MIXED;
        }

        if(startContext == "")
            startContext = CONTEXT_INTRONIC;

        if(endContext == "")
            endContext = CONTEXT_INTRONIC;

        AltSpliceJunctionType type = CRYPTIC_SPLICE_SITE;

        // check for skipped exons indicated by the same transcript matching the start and end of the alt SJ, but skipping an exon
        long nearestStartExon = !sjStartRegions.isEmpty() ? 0 : -1;
        long nearestEndExon = !sjEndRegions.isEmpty() ? 0 : -1;

        if (!sjStartRegions.isEmpty() || !sjEndRegions.isEmpty())
        {
            for (Integer transId : relatedTransIds)
            {
                boolean matchesStart = sjStartRegions.stream().anyMatch(x -> x.hasTransId(transId));
                boolean matchesEnd = sjEndRegions.stream().anyMatch(x -> x.hasTransId(transId));
                if (matchesStart && matchesEnd)
                {
                    type = SKIPPED_EXONS;
                    break;
                }
            }

            if(type != SKIPPED_EXONS)
                type = NEW_EXONS;
        }
        else
        {
            for(RegionReadData region : mGene.getExonRegions())
            {
                if(nearestStartExon != 0)
                {
                    long exonDistance = abs(region.end() - spliceJunction[SE_START]);
                    nearestStartExon = nearestStartExon > 0 ? min(nearestStartExon, exonDistance) : exonDistance;
                }

                if(nearestEndExon != 0)
                {
                    long exonDistance = abs(region.start() - spliceJunction[SE_END]);
                    nearestEndExon = nearestEndExon > 0 ? min(nearestEndExon, exonDistance) : exonDistance;
                }
            }
        }

        AltSpliceJunction altSplicJunction = new AltSpliceJunction(mGene, spliceJunction, type, startContext, endContext);

        altSplicJunction.StartTranscripts.addAll(sjStartRegions);
        altSplicJunction.EndTranscripts.addAll(sjEndRegions);

        return altSplicJunction;
    }

    public void registerAltSpliceJunction(final ReadRecord read, final List<Integer> invalidTranscripts)
    {
        if(mGene == null)
            return;

        AltSpliceJunction altSpliceJunc = createFromRead(read, invalidTranscripts);

        AltSpliceJunction existingSpliceJunc = mAltSpliceJunctions.stream()
                .filter(x -> x.matches(altSpliceJunc)).findFirst().orElse(null);

        if(existingSpliceJunc != null)
        {
            existingSpliceJunc.addFragmentCount();
            return;
        }

        // extra check that the SJ doesn't match any transcript
        if(junctionMatchesGene(altSpliceJunc.SpliceJunction, invalidTranscripts))
            return;

        altSpliceJunc.addFragmentCount();
        mAltSpliceJunctions.add(altSpliceJunc);
    }

    public void recordDepthCounts()
    {
        BamSlicer slicer = new BamSlicer(DEFAULT_MIN_MAPPING_QUALITY, true);

        QueryInterval[] queryInterval = new QueryInterval[1];
        int chrSeqIndex = mSamReader.getFileHeader().getSequenceIndex(mGene.GeneData.Chromosome);

        for(AltSpliceJunction altSJ : mAltSpliceJunctions)
        {
            for(int se = SE_START; se <= SE_END; ++se)
            {
                int position = (int)altSJ.SpliceJunction[se];
                queryInterval[0] = new QueryInterval(chrSeqIndex, position, position);
                int depth = slicer.slice(mSamReader, queryInterval).size();
                altSJ.setFragmentCount(se, depth);
            }
        }
    }

    public static BufferedWriter createAltSpliceJunctionWriter(final RnaExpConfig config)
    {
        try
        {
            final String outputFileName = config.formOutputFile("alt_splice_junc.csv");

            BufferedWriter writer = createBufferedWriter(outputFileName, false);
            writer.write("GeneId,GeneName,Chromosome,Strand,SjStart,SjEnd,FragCount");
            writer.write(",Type,StartContext,EndContext,NearestStartExon,NearestEndExon");
            writer.write(",StartBases,EndBases,StartTrans,EndTrans");
            writer.newLine();
            return writer;
        }
        catch (IOException e)
        {
            RE_LOGGER.error("failed to create at splice junction writer: {}", e.toString());
            return null;
        }
    }

    public void writeAltSpliceJunctions()
    {
        writeAltSpliceJunctions(mWriter, mAltSpliceJunctions, mConfig);
    }

    private synchronized static void writeAltSpliceJunctions(
            final BufferedWriter writer, final List<AltSpliceJunction> altSpliceJunctions, final RnaExpConfig config)
    {
        try
        {
            for(final AltSpliceJunction altSJ : altSpliceJunctions)
            {
                writer.write(String.format("%s,%s,%s,%d",
                        altSJ.Gene.GeneData.GeneId, altSJ.Gene.GeneData.GeneName,
                        altSJ.Gene.GeneData.Chromosome, altSJ.Gene.GeneData.Strand));

                writer.write(String.format(",%d,%d,%d",
                        altSJ.SpliceJunction[SE_START], altSJ.SpliceJunction[SE_END], altSJ.getFragmentCount()));

                final String startTrans = altSJ.StartTranscripts.isEmpty() ? "NONE" : altSJ.startTranscriptNames();
                final String endTrans = altSJ.EndTranscripts.isEmpty() ? "NONE" : altSJ.endTranscriptNames();

                writer.write(String.format(",%s,%s,%s,%d,%d,%s,%s,%s,%s",
                        altSJ.Type, altSJ.RegionContexts[SE_START], altSJ.RegionContexts[SE_END],
                        altSJ.nearestStartExon(), altSJ.nearestEndExon(),
                        altSJ.getBaseContext(config.RefFastaSeqFile, SE_START), altSJ.getBaseContext(config.RefFastaSeqFile, SE_END),
                        startTrans, endTrans));

                writer.newLine();
            }

        }
        catch(IOException e)
        {
            RE_LOGGER.error("failed to write alt splice junction file: {}", e.toString());
        }
    }



}
