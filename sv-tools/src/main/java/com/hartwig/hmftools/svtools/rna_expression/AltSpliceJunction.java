package com.hartwig.hmftools.svtools.rna_expression;

import static java.lang.Math.abs;
import static java.lang.Math.min;

import static com.hartwig.hmftools.linx.analysis.SvUtilities.appendStrList;
import static com.hartwig.hmftools.linx.types.SvVarData.SE_END;
import static com.hartwig.hmftools.linx.types.SvVarData.SE_PAIR;
import static com.hartwig.hmftools.linx.types.SvVarData.SE_START;
import static com.hartwig.hmftools.sig_analyser.loaders.SigSnvLoader.convertBase;
import static com.hartwig.hmftools.svtools.rna_expression.AltSpliceJunctionType.CRYPTIC_SPLICE_SITE;
import static com.hartwig.hmftools.svtools.rna_expression.AltSpliceJunctionType.NEW_EXONS;
import static com.hartwig.hmftools.svtools.rna_expression.AltSpliceJunctionType.SKIPPED_EXONS;
import static com.hartwig.hmftools.svtools.rna_expression.RegionReadData.extractTransName;
import static com.hartwig.hmftools.svtools.rna_expression.RnaExpUtils.positionWithin;
import static com.hartwig.hmftools.svtools.rna_expression.TransMatchType.SPLICE_JUNCTION;

import java.util.List;
import java.util.Map;
import java.util.stream.Collectors;

import com.google.common.collect.Lists;

import htsjdk.samtools.CigarOperator;
import htsjdk.samtools.reference.IndexedFastaSequenceFile;

public class AltSpliceJunction
{
    public final GeneReadData Gene;
    public final long[] SpliceJunction;
    public final AltSpliceJunctionType Type;
    public final List<RegionReadData> StartTranscripts;
    public final List<RegionReadData> EndTranscripts;

    public final String[] RegionContexts;

    public static final String CONTEXT_SJ = "SJ";
    public static final String CONTEXT_EXONIC = "EXONIC";
    public static final String CONTEXT_INTRONIC = "INTRONIC";
    public static final String CONTEXT_MIXED = "MIXED";

    private int mFragmentCount;

    /*
    Record each novel splice junction per Gene + following fields:
        gene, chromosome, pos start & end (ie splice junction)
        count of observations
        list of transcripts that share an exon boundary at start
        list of transcripts that share an exon boundary at end
        distance to nearest known splice boundary at start
        distance to nearest known splice boundary at end
        Start count per (transcript combination) category
        End count per (transcript combination) category
        Annotate skipped exons, cassette exons, cryptic splice sites
        Use to analyse AHR and APC novel splice junctions in relevant samples.
        Generate PON
        Look for retained introns

     */

    public AltSpliceJunction(
            final GeneReadData geneReadData, final long[] spliceJunction, AltSpliceJunctionType type,
            final String startContext, final String endContext)
    {
        Gene = geneReadData;
        SpliceJunction = spliceJunction;

        RegionContexts = new String[] {startContext, endContext};

        mFragmentCount = 0;
        StartTranscripts = Lists.newArrayList();
        EndTranscripts = Lists.newArrayList();

        Type = type;
    }

    public static boolean isCandidate(final ReadRecord read)
    {
        if(!read.Cigar.containsOperator(CigarOperator.N))
            return false;

        if(read.getTranscriptClassifications().values().contains(SPLICE_JUNCTION))
            return false;

        return read.getMappedRegionCoords().size() >= 2;
    }

    public boolean matches(final AltSpliceJunction other)
    {
        if(!other.Gene.GeneData.GeneId.equals(Gene.GeneData.GeneId))
            return false;

        return SpliceJunction[SE_START] == other.SpliceJunction[SE_START] && SpliceJunction[SE_END] == other.SpliceJunction[SE_END];
    }

    public static AltSpliceJunction fromRead(
            final GeneReadData geneReadData, final ReadRecord read, final List<Integer> relatedTransIds)
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
            for(RegionReadData region : geneReadData.getExonRegions())
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

        AltSpliceJunction altSplicJunction = new AltSpliceJunction(geneReadData, spliceJunction, type, startContext, endContext);

        altSplicJunction.StartTranscripts.addAll(sjStartRegions);
        altSplicJunction.EndTranscripts.addAll(sjEndRegions);

        return altSplicJunction;
    }

    public int getFragmentCount() { return mFragmentCount;}
    public void addFragmentCount() { ++mFragmentCount;}

    public String startTranscriptNames() { return generateTranscriptNames( StartTranscripts); }
    public String endTranscriptNames() { return generateTranscriptNames( EndTranscripts); }

    private String generateTranscriptNames(final List<RegionReadData> regions)
    {
        List<String> transNames = Lists.newArrayList();

        for(RegionReadData region: regions)
        {
            transNames.addAll(region.getRefRegions().stream().map(x -> extractTransName(x)).collect(Collectors.toList()));
        }

        return appendStrList(transNames, ';');
    }

    public long nearestStartExon()
    {
        if(!StartTranscripts.isEmpty())
            return 0;

        long nearestStartExon = -1;

        for(RegionReadData region : Gene.getExonRegions())
        {
            if(nearestStartExon != 0)
            {
                long exonDistance = abs(region.end() - SpliceJunction[SE_START]);
                nearestStartExon = nearestStartExon > 0 ? min(nearestStartExon, exonDistance) : exonDistance;
            }
        }

        return nearestStartExon;
    }

    public long nearestEndExon()
    {
        if(!EndTranscripts.isEmpty())
            return 0;

        long nearestEndExon = -1;

        for(RegionReadData region : Gene.getExonRegions())
        {
            if(nearestEndExon != 0)
            {
                long exonDistance = abs(region.start() - SpliceJunction[SE_END]);
                nearestEndExon = nearestEndExon > 0 ? min(nearestEndExon, exonDistance) : exonDistance;
            }
        }

        return nearestEndExon;
    }

    public String getBaseContext(final IndexedFastaSequenceFile refGenome, int seIndex)
    {
        long position = SpliceJunction[seIndex];
        int startOffset = ((seIndex == SE_START) == (Gene.GeneData.Strand == 1)) ? 1 : 10;
        int endOffset = startOffset == 1 ? 10: 1;

        final String baseStr = refGenome
                .getSubsequenceAt(Gene.GeneData.Chromosome, position - startOffset, position + endOffset).getBaseString();

        if(Gene.GeneData.Strand == 1)
            return baseStr;

        String reverseStr = "";

        for(int i = baseStr.length() - 1; i >= 0; --i)
        {
            reverseStr += String.valueOf(convertBase(baseStr.charAt(i)));
        }

        return reverseStr;
    }



}
